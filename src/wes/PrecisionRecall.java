package wes;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import utils.Region;

public class PrecisionRecall {

    private static final Map<String, Set<Region>> knownCNVMap         = new HashMap<>();
    private static final Map<String, Set<Region>> seqcnvPredictMap    = new HashMap<>();
    private static final Map<String, Set<Region>> coniferPredictMap   = new HashMap<>();
    private static final Map<String, Set<Region>> cnvnatorPredictMap  = new HashMap<>();
    private static final Map<String, Set<Region>> xhmmPredictMap      = new HashMap<>();
    private static final Map<String, Set<Region>> excavatorPredictMap = new HashMap<>();
    private static final Set<String>              excludedSamples     = new HashSet<>(
        Arrays.asList("NA18959", "NA18999", "NA19152"));

    public static void main(String[] args) {

        System.out.println("WES data precision recall analysis:");
        System.out.println("Recall = correctly detected events / known CNV events");
        System.out.println("Precision = unique correctly detected events / total tool CNV events");
        System.out.println();

        float[] overlapRatios = { 0.1f };
        for (float overlapRatio : overlapRatios) {
            seqcnv(overlapRatio, "C:\\Users\\Administrator\\Desktop\\chr1_chr4_result\\SeqCNV",
                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv");
            conifer(overlapRatio,
                "C:\\Users\\Administrator\\Desktop\\chr1_chr4_result\\CoNIFER\\CoNIFER_chr1_chr4.tsv",
                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv");
            cnvnator(overlapRatio, "C:\\Users\\Administrator\\Desktop\\chr1_chr4_result\\CNVnator",
                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv");
            xhmm(overlapRatio,
                "C:\\Users\\Administrator\\Desktop\\chr1_chr4_result\\XHMM\\DATA.xcnv",
                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv");
            excavator(overlapRatio,
                "C:\\Users\\Administrator\\Desktop\\chr1_chr4_result\\EXCAVATOR",
                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv");
        }
    }

    /**
     * SeqCNV的overlap分析，输出precision和recall。
     * <pre>
     * Note:
     * 1. 没有考虑SeqCNV的CNV region可能会包含多个known CNV region的情况。
     * 2. 如果known CNV region overlap了好几个predict CNV region, 会导致重复计算。
     * 
     * Problem：
     * knownCNV中重复的region要不要去掉？
     * 1. 不去掉，干脆重复计算。
     * 2. 去掉更短的重复region。（待做）
     * </pre>
     * 
     * @param overlapRatio
     * @param seqcnvResultFolder
     * @param knownCNVFolder
     */
    public static void seqcnv(float overlapRatio, String seqcnvResultFolder,
                              String knownCNVFolder) {
        System.out.println("SeqCNV:");
        System.out.println("Overlap ratio threshold is " + overlapRatio);

        File folderFile = new File(seqcnvResultFolder);
        String[] sampleNames = folderFile.list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                if (name.startsWith("NA"))
                    return true;
                return false;
            }
        });

        for (String sample : sampleNames) {
            if (excludedSamples.contains(sample))
                continue; // 如果在排除样本里
            String resultFilePath = seqcnvResultFolder + File.separator + sample + File.separator
                                    + "report" + File.separator + "CNV_report.txt";
            Set<Region> seqcnvRegions = new HashSet<>();
            List<String> seqcnvResultLines = readLines(resultFilePath);
            seqcnvResultLines.stream().forEach(line -> {
                String[] feature = line.split("\\s+");
                seqcnvRegions
                    .add(new Region(feature[0].length() < 3 ? "chr" + feature[0] : feature[0],
                        feature[1], feature[2]));
            });
            seqcnvPredictMap.put(sample, seqcnvRegions);
        }
        outputPrecisionRecall(overlapRatio, seqcnvPredictMap, knownCNVFolder);
        System.out.println();
    }

    /**
     * CoNIFER的overlap分析，输出precision和recall。
     * <pre>
     * Note:
     * 1. 没有考虑CoNIFER的CNV region可能会包含多个known CNV region的情况。
     * 2. 如果known CNV region overlap了好几个predict CNV region, 会导致重复计算。
     * 
     * Problem：
     * knownCNV中重复的region要不要去掉？
     * 1. 不去掉，干脆重复计算。
     * 2. 去掉更短的重复region。（待做）
     * </pre>
     * @param overlapRatio
     * @param coniferResultFilePath
     * @param knownCNVFolder
     */
    public static void conifer(float overlapRatio, String coniferResultFilePath,
                               String knownCNVFolder) {
        System.out.println("CoNIFER:");
        System.out.println("Overlap ratio threshold is " + overlapRatio);

        List<String> coniferResultLines = readLines(coniferResultFilePath);
        coniferResultLines.stream().forEach(line -> {
            String[] feature = line.split("\\s+");
            String coniferSample = feature[0];
            if (excludedSamples.contains(coniferSample))
                return; // 如果在排除样本里
            Region region = new Region(feature[1].length() < 3 ? "chr" + feature[1] : feature[1],
                feature[2], feature[3]);
            if (coniferPredictMap.containsKey(coniferSample)) {
                coniferPredictMap.get(coniferSample).add(region);
            } else {
                Set<Region> regions = new HashSet<>();
                regions.add(region);
                coniferPredictMap.put(coniferSample, regions);
            }
        });

        outputPrecisionRecall(overlapRatio, coniferPredictMap, knownCNVFolder);
        System.out.println();
    }

    /**
     * CNVnator的overlap分析，输出precision和recall。
     * <pre>
     * Note:
     * 1. 没有考虑CNVnator的CNV region可能会包含多个known CNV region的情况。
     * 2. 如果known CNV region overlap了好几个predict CNV region, 会导致重复计算。
     * 
     * Problem：
     * knownCNV中重复的region要不要去掉？
     * 1. 不去掉，干脆重复计算。
     * 2. 去掉更短的重复region。（待做）
     * </pre>
     * @param overlapRatio
     * @param cnvnatorResultFolder
     * @param knownCNVFolder
     */
    public static void cnvnator(float overlapRatio, String cnvnatorResultFolder,
                                String knownCNVFolder) {
        System.out.println("CNVnator:");
        System.out.println("Overlap ratio threshold is " + overlapRatio);

        File folderFile = new File(cnvnatorResultFolder);
        String[] sampleNames = folderFile.list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                if (name.startsWith("NA"))
                    return true;
                return false;
            }
        });
        for (String sample : sampleNames) {
            if (excludedSamples.contains(sample))
                continue; // 如果在排除样本里
            String resultFilePath = cnvnatorResultFolder + File.separator + sample + File.separator
                                    + sample + ".call";
            List<String> cnvnatorResultLines = readLines(resultFilePath);
            Set<Region> cnvnatorRegions = new HashSet<>();
            cnvnatorResultLines.stream().forEach(line -> {
                String[] feature = line.split("\\s+")[1].replace(":", "-").split("-");
                cnvnatorRegions
                    .add(new Region(feature[0].length() < 3 ? "chr" + feature[0] : feature[0],
                        feature[1], feature[2]));
            });
            cnvnatorPredictMap.put(sample, cnvnatorRegions);
        }

        outputPrecisionRecall(overlapRatio, cnvnatorPredictMap, knownCNVFolder);
        System.out.println();
    }

    /**
     * XHMM的overlap分析，输出precision和recall。
     * <pre>
     * Note:
     * 1. 没有考虑XHMM的CNV region可能会包含多个known CNV region的情况。
     * 2. 如果known CNV region overlap了好几个predict CNV region, 会导致重复计算。
     * 
     * Problem：
     * knownCNV中重复的region要不要去掉？
     * 1. 不去掉，干脆重复计算。
     * 2. 去掉更短的重复region。（待做）
     * </pre>
     * 
     * @param overlapRatio
     * @param xhmmResultFilePath
     * @param knownCNVFolder
     */
    public static void xhmm(float overlapRatio, String xhmmResultFilePath, String knownCNVFolder) {
        System.out.println("XHMM:");
        System.out.println("Overlap ratio threshold is " + overlapRatio);

        List<String> xhmmResultLines = readLines(xhmmResultFilePath);
        xhmmResultLines.remove(0); // drop header
        xhmmResultLines.stream().forEach(line -> {
            String[] feature = line.split("\\s+");
            String sample = feature[0];
            if (excludedSamples.contains(sample))
                return; // 如果在排除样本里
            String[] regionFeature = feature[2].replace(":", "-").split("-");
            Region region = new Region(
                regionFeature[0].length() < 3 ? "chr" + regionFeature[0] : regionFeature[0],
                regionFeature[1], regionFeature[2]);
            if (xhmmPredictMap.containsKey(sample)) {
                xhmmPredictMap.get(sample).add(region);
            } else {
                Set<Region> regions = new HashSet<>();
                regions.add(region);
                xhmmPredictMap.put(sample, regions);
            }
        });

        outputPrecisionRecall(overlapRatio, xhmmPredictMap, knownCNVFolder);
        System.out.println();
    }

    /**
     * EXCAVATOR的overlap分析，输出precision和recall。
     * <pre>
     * Note:
     * 1. 没有考虑EXCAVATOR的CNV region可能会包含多个known CNV region的情况。
     * 2. 如果known CNV region overlap了好几个predict CNV region, 会导致重复计算。
     * 
     * Problem：
     * knownCNV中重复的region要不要去掉？
     * 1. 不去掉，干脆重复计算。
     * 2. 去掉更短的重复region。（待做）
     * </pre>
     * 
     * @param overlapRatio
     * @param excavatorResultFolder
     * @param knownCNVFolder
     */
    public static void excavator(float overlapRatio, String excavatorResultFolder,
                                 String knownCNVFolder) {
        System.out.println("EXCAVATOR:");
        System.out.println("Overlap ratio threshold is " + overlapRatio);

        File folderFile = new File(excavatorResultFolder);
        String[] sampleNames = folderFile.list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                if (name.startsWith("NA"))
                    return true;
                return false;
            }
        });

        for (String sample : sampleNames) {
            if (excludedSamples.contains(sample))
                continue; // 如果在排除样本里
            String resultFilePath = excavatorResultFolder + File.separator + sample + File.separator
                                    + "FastCallResults_" + sample + ".txt";
            List<String> excavatorResultLines = readLines(resultFilePath);
            excavatorResultLines.remove(0);
            Set<Region> excavatorRegions = new HashSet<>();
            excavatorResultLines.stream().forEach(line -> {
                String[] feature = line.split("\\s+");
                excavatorRegions
                    .add(new Region(feature[0].length() < 3 ? "chr" + feature[0] : feature[0],
                        feature[1], feature[2]));
            });
            excavatorPredictMap.put(sample, excavatorRegions);
        }

        outputPrecisionRecall(overlapRatio, excavatorPredictMap, knownCNVFolder);
        System.out.println();
    }

    /**
     * read known CNV regions and calculate precision & recall.
     * @param overlapRatio
     * @param toolPredictMap
     * @param knownCNVFolder
     */
    public static void outputPrecisionRecall(Float overlapRatio,
                                             Map<String, Set<Region>> toolPredictMap,
                                             String knownCNVFolder) {

        Set<Region> uniquePredictedRegions = new HashSet<>();

        int correctedCNVRegions = 0;
        int uniqueCorrectedCNVRegions = 0;
        int knownCNVRegions = 0;
        int toolCNVRegions = 0;

        File knownCNVFolderFile = new File(knownCNVFolder);
        String[] knownCNVFileNames = knownCNVFolderFile.list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                if (name.endsWith("tsv"))
                    return true;
                return false;
            }
        });

        for (String knownCNVFileName : knownCNVFileNames) {
            String[] tmps = knownCNVFileName.replace("_", ".").split("\\."); // split里默认是正则表达式，要按照.分割的话需要转义
            String sample = tmps[tmps.length - 2];
            if (excludedSamples.contains(sample))
                continue;
            String knownCNVFilePath = knownCNVFolder + File.separator + knownCNVFileName;
            readKnownCNVRegions(knownCNVFilePath, sample);
        }

        toolCNVRegions += toolPredictMap.entrySet().stream().mapToInt(entry -> {
            return entry.getValue().size();
        }).sum();

        knownCNVRegions += knownCNVMap.entrySet().stream().mapToInt(entry -> {
            String sample = entry.getKey();
            if (excludedSamples.contains(sample))
                return 0;
            Set<Region> knownRegions = entry.getValue();
            return knownRegions.stream().mapToInt(region -> {
                if (region.getChr().equals("chr1") || region.getChr().equals("chr4"))
                    return 1;
                return 0;
            }).sum();
        }).sum();

        for (Map.Entry<String, Set<Region>> entry : toolPredictMap.entrySet()) {
            Set<Region> seqRegions = entry.getValue();
            Set<Region> knownRegions = knownCNVMap.get(entry.getKey());

            int predictCount = 0;
            for (Region predictRegion : seqRegions) {
                for (Region knownRegion : knownRegions) {
                    if (knownRegion.intersact(predictRegion)) {
                        long predictLength = knownRegion.overlap(predictRegion);
                        if ((float) predictLength
                            / (float) knownRegion.getLength() >= overlapRatio) {
                            predictCount++;
                            uniquePredictedRegions.add(predictRegion);
                        }
                    }
                }
            }
            correctedCNVRegions += predictCount;
        }
        uniqueCorrectedCNVRegions = uniquePredictedRegions.size();

        System.out.println("Correctly detected CNV events: " + correctedCNVRegions
                           + ", Unique correctly detected CNV events: " + uniqueCorrectedCNVRegions
                           + ", Known CNV events: " + knownCNVRegions + ", Tool CNV events: "
                           + toolCNVRegions);
        System.out
            .println(
                "Recall: " + String.format("%.2f", (double) correctedCNVRegions / knownCNVRegions)
                     + ", Precision: "
                     + String.format("%.2f", (double) uniqueCorrectedCNVRegions / toolCNVRegions));
    }

    /**
     * sum the overlap between the knownCNVRegion and predictCNVRegion, return the num of
     * knownCNVRegions which was overlapped by predictCNVRegions exceed the overlapRatio.
     * 
     * @param beneathRatioMap key is knownCNVRegion, value is a set of predictCNVRegions which
     *        overlap with the key but beneath the overlapRatio.
     * @param overlapRatio
     * @return num of the overlap regions exceed the overlapRatio
     */
    private static int sumOverlap(Map<Region, Set<Region>> beneathRatioMap, float overlapRatio) {
        int sum = 0;
        for (Map.Entry<Region, Set<Region>> entry : beneathRatioMap.entrySet()) {
            Region knownRegion = entry.getKey();
            Set<Region> predictRegions = entry.getValue();
            long overlapLength = 0;
            for (Region predictRegion : predictRegions) {
                overlapLength += knownRegion.overlap(predictRegion);
            }
            if ((double) overlapLength / knownRegion.getLength() >= overlapRatio) {
                sum++;
            }
        }
        return sum;
    }

    private static void readKnownCNVRegions(String knownCNVFilePath, String sample) {
        if (knownCNVMap.containsKey(sample)) {
            return;
        }
        List<String> knownCNVLines = readLines(knownCNVFilePath);
        Set<Region> knownCNVRegions = new LinkedHashSet<>();
        knownCNVLines.remove(0);
        knownCNVLines.stream().forEach(line -> {
            String[] feature = line.split("\\s+");
            knownCNVRegions.add(new Region(
                feature[0].length() < 3 ? "chr" + feature[0] : feature[0], feature[1], feature[2]));
        });
        knownCNVMap.put(sample, knownCNVRegions);
    }

    /**
     * read lines from file.
     * 
     * @param filePath
     * @return lines of the file.
     */
    private static List<String> readLines(String filePath) {
        List<String> lines = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(filePath))) {
            String line = null;
            while ((line = br.readLine()) != null) {
                lines.add(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return lines;
    }
}