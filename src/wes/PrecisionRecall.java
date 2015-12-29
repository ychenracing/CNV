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

import utils.Pair;
import utils.Region;

/**
 * 
 * @author Administrator
 * @version $Id: PrecisionRecall.java, v 0.1 2015年12月18日 下午1:32:31 Administrator Exp $
 */
public class PrecisionRecall {

    private static final Map<String, Set<Region>> knownCNVMap         = new HashMap<>();
    private static final Map<String, Set<Region>> seqcnvPredictMap    = new HashMap<>();
    private static final Map<String, Set<Region>> coniferPredictMap   = new HashMap<>();
    private static final Map<String, Set<Region>> cnvnatorPredictMap  = new HashMap<>();
    private static final Map<String, Set<Region>> xhmmPredictMap      = new HashMap<>();
    private static final Map<String, Set<Region>> excavatorPredictMap = new HashMap<>();
    private static final Set<String>              sampleSet           = new HashSet<>(
        Arrays.asList("NA10847", "NA12760", "NA11840", "NA12761", "NA12249", "NA18959", "NA12717",
            "NA18966", "NA12751", "NA18967", "NA18970", "NA19138", "NA18973", "NA19153", "NA18981",
            "NA19159", "NA18999", "NA19206", "NA19131", "NA19223"));
    //    private static final Set<String>              excludedSamples     = new HashSet<>(
    //        Arrays.asList("NA18959", "NA18999", "NA19152", "NA18973", "NA12760", "NA18966", "NA19223")); // 前3个影响结果，后几个还未完成
    private static final Set<String>              excludedSamples     = new HashSet<>(
        Arrays.asList("NA19152"));                                                                 // 前3个影响结果，后几个还未完成
    private static final Set<String>              analysedChromosomes = new HashSet<>(
        Arrays.asList("chr1", "chr4"));

    public static void main(String[] args) {

        System.out.println("WES data precision recall analysis:");
        System.out.println("Recall = predicted known CNV events / total known CNV events");
        System.out.println("Precision = unique correctly detected events / total tool CNV events");
        System.out.println("Analysed chromosomes are " + analysedChromosomes);
        System.out.println();

        float[] overlapRatios = { 0.1f };
        for (float overlapRatio : overlapRatios) {
            seqcnv(overlapRatio, "C:\\Users\\Administrator\\Desktop\\chr1_chr4_result\\SeqCNV",
                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv"); // chr1_chr4
            //            seqcnv(overlapRatio, "C:\\Users\\Administrator\\Desktop\\seqcnv_chr1",
            //                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv"); // only chr1
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
            System.out.println();
            System.out.println();
            System.out.println();
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
                feature[0] = feature[0].length() < 3 ? "chr" + feature[0] : feature[0];
                if (analysedChromosomes.contains(feature[0]))
                    seqcnvRegions.add(new Region(feature[0], feature[1], feature[2]));
            });
            seqcnvPredictMap.put(sample, seqcnvRegions);
        }
        //        for (Map.Entry<String, Set<Region>> entry : seqcnvPredictMap.entrySet()) {
        //            System.out.println(entry.getKey() + ":" + entry.getValue().size());
        //        }
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
            feature[1] = feature[1].length() < 3 ? "chr" + feature[1] : feature[1];
            if (!analysedChromosomes.contains(feature[1]))
                return;
            Region region = new Region(feature[1], feature[2], feature[3]);
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
                feature[0] = feature[0].length() < 3 ? "chr" + feature[0] : feature[0];
                if (analysedChromosomes.contains(feature[0]))
                    cnvnatorRegions.add(new Region(feature[0], feature[1], feature[2]));
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
            regionFeature[0] = regionFeature[0].length() < 3 ? "chr" + regionFeature[0]
                : regionFeature[0];
            if (!analysedChromosomes.contains(regionFeature[0]))
                return;
            Region region = new Region(regionFeature[0], regionFeature[1], regionFeature[2]);
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
                feature[0] = feature[0].length() < 3 ? "chr" + feature[0] : feature[0];
                if (analysedChromosomes.contains(feature[0]))
                    excavatorRegions.add(new Region(feature[0], feature[1], feature[2]));
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
        int predictedKnownRegions = 0;
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

        Map<String, Map<Region, Set<Region>>> beneathMap = new HashMap<>(); // key为sample
        Map<String, Set<Region>> predictedKnownRegionMap = new HashMap<>(); // key为sample
        Map<String, Set<Region>> uniquePredictedRegionMap = new HashMap<>(); // key为sample

        Map<String, Pair<Float, Float>> preRecMap = new HashMap<>(); // key为sample，value为该sample对应的recall和precision。

        for (String sample : sampleSet) {
            if (excludedSamples.contains(sample) || !toolPredictMap.containsKey(sample))
                continue;
            Set<Region> sampleToolCNVRegions = toolPredictMap.get(sample);
            if (sampleToolCNVRegions == null || sampleToolCNVRegions.isEmpty())
                continue;
            Set<Region> sampleKnownCNVRegions = knownCNVMap.get(sample);
            Set<Region> sampleCorrectedCNVRegions = new HashSet<>();
            Set<Region> samplePredictedKnownRegions = new HashSet<>();

            for (Region predictRegion : sampleToolCNVRegions) {
                for (Region knownRegion : sampleKnownCNVRegions) {
                    if (knownRegion.isOverlapped(predictRegion)) {
                        long predictLength = knownRegion.getOverlapLength(predictRegion);

                        //                        if ((float) predictLength
                        //                            / (float) knownRegion.getLength() >= overlapRatio) {

                        if ((double) predictLength / (double) predictRegion
                            .getOverlapBaseLength(knownRegion) >= overlapRatio) {
                            //                            System.out.println(predictLength + ":"
                            //                                               + predictRegion.getOverlapBaseLength(knownRegion));
                            //                            if (predictRegion.getOverlapBaseLength(knownRegion) == 0) {
                            //                                System.out.println(predictRegion + " " + knownRegion);
                            //                            }
                            samplePredictedKnownRegions.add(knownRegion);
                            sampleCorrectedCNVRegions.add(predictRegion);

                            if (uniquePredictedRegionMap.containsKey(sample)) {
                                uniquePredictedRegionMap.get(sample).add(predictRegion);
                            } else {
                                Set<Region> predicted = new HashSet<>();
                                predicted.add(predictRegion);
                                uniquePredictedRegionMap.put(sample, predicted);
                            }
                            if (predictedKnownRegionMap.containsKey(sample)) {
                                predictedKnownRegionMap.get(sample).add(knownRegion);
                            } else {
                                Set<Region> predicted = new HashSet<>();
                                predicted.add(knownRegion);
                                predictedKnownRegionMap.put(sample, predicted);
                            }
                        }
                    } else {
                        if (beneathMap.containsKey(sample)) {
                            Map<Region, Set<Region>> beneathMapItem = beneathMap.get(sample);
                            if (beneathMapItem.containsKey(knownRegion)) {
                                beneathMapItem.get(knownRegion).add(predictRegion);
                            } else {
                                Set<Region> beneathSet = new HashSet<>();
                                beneathSet.add(predictRegion);
                                beneathMapItem.put(knownRegion, beneathSet);
                            }
                        } else {
                            Set<Region> beneathSet = new HashSet<>();
                            beneathSet.add(predictRegion);
                            Map<Region, Set<Region>> beneathMapItem = new HashMap<>();
                            beneathMapItem.put(knownRegion, beneathSet);
                            beneathMap.put(sample, beneathMapItem);
                        }
                    }
                }
            }

            int extraPredictedKnown = 0;
            Map<Region, Set<Region>> sampleBeneathMap = beneathMap.get(sample);
            if (sampleBeneathMap != null && !sampleBeneathMap.isEmpty()) {
                for (Map.Entry<Region, Set<Region>> beneathEntry : sampleBeneathMap.entrySet()) {
                    Region knownRegion = beneathEntry.getKey();
                    Set<Region> predictRegions = beneathEntry.getValue();
                    long overlapLength = 0;
                    for (Region predictRegion : predictRegions) {
                        overlapLength += knownRegion.getOverlapLength(predictRegion);
                    }
                    if ((double) overlapLength / knownRegion.getLength() >= overlapRatio) {
                        extraPredictedKnown++;
                    }
                }
            }

            float sampleRecall = (float) (samplePredictedKnownRegions.size() + extraPredictedKnown)
                                 / sampleKnownCNVRegions.size();
            float samplePrecision = (float) sampleCorrectedCNVRegions.size()
                                    / sampleToolCNVRegions.size();
            System.out.println(sample + ":" + String.format("%.2f", sampleRecall) + ", "
                               + String.format("%.2f", samplePrecision));
            preRecMap.put(sample, new Pair<>(sampleRecall, samplePrecision));

        }

        predictedKnownRegions += predictedKnownRegionMap.entrySet().stream().mapToInt(entry ->

        {
            return entry.getValue().size();
        }).sum();

        for (

        Map.Entry<String, Map<Region, Set<Region>>> entry : beneathMap.entrySet())

        {
            Map<Region, Set<Region>> beneathMapItem = entry.getValue();
            for (Map.Entry<Region, Set<Region>> beneathEntry : beneathMapItem.entrySet()) {
                Region knownRegion = beneathEntry.getKey();
                Set<Region> predictRegions = beneathEntry.getValue();
                long overlapLength = 0;
                for (Region predictRegion : predictRegions) {
                    overlapLength += knownRegion.getOverlapLength(predictRegion);
                }
                if ((double) overlapLength / knownRegion.getLength() >= overlapRatio) {
                    predictedKnownRegions++;
                }
            }
        }

        uniqueCorrectedCNVRegions = uniquePredictedRegionMap.entrySet().stream().mapToInt(entry ->

        {
            return entry.getValue().size();
        }).sum();

        System.out.println("Predicted known CNV events: " + predictedKnownRegions
                           + ", Unique correctly detected CNV events: " + uniqueCorrectedCNVRegions
                           + ", Known CNV events: " + knownCNVRegions + ", Tool CNV events: "
                           + toolCNVRegions);
        System.out
            .println(
                "Recall: " + String.format("%.2f", (double) predictedKnownRegions / knownCNVRegions)
                     + ", Precision: "
                     + String.format("%.2f", (double) uniqueCorrectedCNVRegions / toolCNVRegions));
        sampleSet.addAll(excludedSamples);
        sampleSet.removeAll(excludedSamples);

        int analysedSampleCount = sampleSet.size();
        double avgRecall = preRecMap.entrySet().stream().mapToDouble(entry -> {
            return entry.getValue().getFirst();
        }).sum() / analysedSampleCount;
        double avgPrecision = preRecMap.entrySet().stream().mapToDouble(entry -> {
            return entry.getValue().getSecond();
        }).sum() / analysedSampleCount;
        System.out.println("Recall: " + String.format("%.2f", avgRecall) + ", Precision: "
                           + String.format("%.2f", avgPrecision) + " (average)");

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
                overlapLength += knownRegion.getOverlapLength(predictRegion);
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
            feature[0] = feature[0].length() < 3 ? "chr" + feature[0] : feature[0];
            if (analysedChromosomes.contains(feature[0]))
                knownCNVRegions.add(new Region(feature[0], feature[1], feature[2]));
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
