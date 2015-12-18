package wes;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import utils.Region;

public class OverlapAnalysis {

    private static final Map<String, Set<Region>> knownCNVMap        = new HashMap<>();
    private static final Map<String, Set<Region>> seqcnvPredictMap   = new HashMap<>();
    private static final Map<String, Set<Region>> coniferPredictMap  = new HashMap<>();
    private static final Map<String, Set<Region>> cnvnatorPredictMap = new HashMap<>();
    private static final Map<String, Set<Region>> xhmmPredictMap     = new HashMap<>();

    public static void main(String[] args) {
        float overlapRatio = 0.1f;
        seqcnv(overlapRatio, "C:\\Users\\Administrator\\Desktop\\chr1_chr4_result\\SeqCNV",
            "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv");
        conifer(overlapRatio,
            "C:\\Users\\Administrator\\Desktop\\chr1_chr4_result\\CoNIFER\\CoNIFER_chr1_chr4.tsv",
            "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv");
        cnvnator(overlapRatio, "C:\\Users\\Administrator\\Desktop\\chr1_chr4_result\\CNVnator",
            "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv");
        xhmm(overlapRatio, "C:\\Users\\Administrator\\Desktop\\chr1_chr4_result\\XHMM\\DATA.xcnv",
            "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv");
    }

    /**
     * SeqCNV的overlap分析。考虑了SeqCNV的CNV region可能会包含多个known CNV region的情况。
     * knownCNV中重复的region要不要去掉？
     * 1. 不去掉，干脆重复计算。
     * 2. 去掉更短的重复region。（待做）
     * 
     * @param overlapRatio
     * @param seqcnvResultFolder
     * @param knownCNVFolder
     */
    public static void seqcnv(float overlapRatio, String seqcnvResultFolder,
                              String knownCNVFolder) {
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
            String resultFilePath = seqcnvResultFolder + File.separator + sample + File.separator
                                    + "report" + File.separator + "CNV_report.txt";
            String knownCNVFilePath = knownCNVFolder + File.separator + "dgv_known_cnv_" + sample
                                      + ".tsv";
            List<String> seqcnvResultLines = readLines(resultFilePath);
            Set<Region> seqcnvRegions = new HashSet<>();
            seqcnvResultLines.stream().forEach(line -> {
                String[] feature = line.split("\\s+");
                seqcnvRegions
                    .add(new Region(feature[0].length() < 3 ? "chr" + feature[0] : feature[0],
                        feature[1], feature[2]));
            });
            seqcnvPredictMap.put(sample, seqcnvRegions);
            readKnownCNVRegions(knownCNVFilePath, sample);
        }

        for (Map.Entry<String, Set<Region>> entry : seqcnvPredictMap.entrySet()) {
            Set<Region> seqRegions = entry.getValue();
            Set<Region> knownRegions = knownCNVMap.get(entry.getKey());
            int predictCount = 0;
            int knownCount = knownRegions.stream().mapToInt(region -> {
                if (region.getChr().equals("chr1") || region.getChr().equals("chr4"))
                    return 1;
                else
                    return 0;
            }).sum();

            // overlap没有超过ratio的CNV region映射，key为knownCNVRegion，value为与该knownCNVRegion overlap的那些predictRegions。
            Map<Region, Set<Region>> beneathRatioMap = new HashMap<>();
            for (Region predictRegion : seqRegions) {

                //                boolean exceedRatio = false;
                for (Region knownRegion : knownRegions) {
                    if (knownRegion.intersact(predictRegion)) {
                        long predictLength = knownRegion.overlap(predictRegion);
                        if ((float) predictLength
                            / (float) knownRegion.getLength() >= overlapRatio) {
                            predictCount++;
                            //                            exceedRatio = true;
                            // System.out.println(knownRegion + " " + predictRegion);
                        } else {
                            if (!beneathRatioMap.containsKey(knownRegion)) {
                                Set<Region> set = new HashSet<>();
                                set.add(predictRegion);
                                beneathRatioMap.put(knownRegion, set);
                            } else {
                                beneathRatioMap.get(knownRegion).add(predictRegion);
                            }
                        }
                    }
                }
                //                if (exceedRatio && beneathRatioMap.containsKey(knownRegion)) {
                //                    beneathRatioMap.remove(knownRegion);
                //                }
            }
            //            predictCount += sumOverlap(beneathRatioMap, overlapRatio);
            System.out.println("SeqCNV: " + entry.getKey() + ": " + predictCount
                               + ", seqcnvPredictRegions: " + seqRegions.size());
        }
    }

    /**
     * CoNIFER的overlap分析。没有考虑CoNIFER的CNV region可能会包含多个known CNV region的情况。
     * 
     * @param overlapRatio
     */
    public static void conifer(float overlapRatio, String coniferResultFilePath,
                               String knownCNVFolder) {
        File knownCNVFolderFile = new File(knownCNVFolder);
        String[] knownCNVFileNames = knownCNVFolderFile.list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                if (name.endsWith("tsv"))
                    return true;
                return false;
            }
        });

        List<String> coniferResultLines = readLines(coniferResultFilePath);
        coniferResultLines.stream().forEach(line -> {
            String[] feature = line.split("\\s+");
            Region region = new Region(feature[1].length() < 3 ? "chr" + feature[1] : feature[1],
                feature[2], feature[3]);
            String coniferSample = feature[0];
            if (coniferPredictMap.containsKey(coniferSample)) {
                coniferPredictMap.get(coniferSample).add(region);
            } else {
                Set<Region> regions = new HashSet<>();
                regions.add(region);
                coniferPredictMap.put(coniferSample, regions);
            }
        });

        for (String knownCNVFileName : knownCNVFileNames) {
            String[] tmps = knownCNVFileName.replace("_", ".").split("\\."); // split里默认是正则表达式，要按照.分割的话需要转义
            String sample = tmps[tmps.length - 2];
            String knownCNVFilePath = knownCNVFolder + File.separator + knownCNVFileName;
            readKnownCNVRegions(knownCNVFilePath, sample);
        }

        for (Map.Entry<String, Set<Region>> entry : coniferPredictMap.entrySet()) {
            if (!knownCNVMap.containsKey(entry.getKey())) {
                continue;
            }
            Set<Region> conRegions = entry.getValue();
            Set<Region> knownRegions = knownCNVMap.get(entry.getKey());
            int predictCount = 0;
            int knownCount = knownRegions.stream().mapToInt(region -> {
                if (region.getChr().equals("chr1") || region.getChr().equals("chr4"))
                    return 1;
                else
                    return 0;
            }).sum();
            Map<Region, Set<Region>> beneathRatioMap = new HashMap<>();
            for (Region predictRegion : conRegions) {
                //                boolean exceedRatio = false;
                for (Region knownRegion : knownRegions) {
                    if (!predictRegion.getChr().equals("chr1")
                        && !predictRegion.getChr().equals("chr4")) {
                        continue; // 只考虑chr1和chr4
                    }
                    if (knownRegion.intersact(predictRegion)) {
                        long predictLength = knownRegion.overlap(predictRegion);
                        if ((double) predictLength
                            / (double) knownRegion.getLength() >= overlapRatio) {
                            predictCount++;
                            //                            exceedRatio = true;
                        } else {
                            if (!beneathRatioMap.containsKey(knownRegion)) {
                                Set<Region> set = new HashSet<>();
                                set.add(predictRegion);
                                beneathRatioMap.put(knownRegion, set);
                            } else {
                                beneathRatioMap.get(knownRegion).add(predictRegion);
                            }
                        }
                    }
                }
                //                if (exceedRatio && beneathRatioMap.containsKey(knownRegion)) {
                //                    beneathRatioMap.remove(knownRegion);
                //                }
            }
            //            predictCount += sumOverlap(beneathRatioMap, overlapRatio);
            if (knownCount != 0) {
                //                System.out.println("CoNIFER: " + entry.getKey() + ": " + predictCount
                //                        + ", known CNV regions: " + knownCount);
                System.out.println("CoNIFER: " + entry.getKey() + ": " + predictCount
                                   + ", coniferPredictRegions: " + conRegions.size());
            }
        }
    }

    public static void cnvnator(float overlapRatio, String cnvnatorResultFolder,
                                String knownCNVFolder) {

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
            String resultFilePath = cnvnatorResultFolder + File.separator + sample + File.separator
                                    + sample + ".call";
            String knownCNVFilePath = knownCNVFolder + File.separator + "dgv_known_cnv_" + sample
                                      + ".tsv";
            List<String> cnvnatorResultLines = readLines(resultFilePath);
            Set<Region> cnvnatorRegions = new HashSet<>();
            cnvnatorResultLines.stream().forEach(line -> {
                String[] feature = line.split("\\s+")[1].replace(":", "-").split("-");
                cnvnatorRegions
                    .add(new Region(feature[0].length() < 3 ? "chr" + feature[0] : feature[0],
                        feature[1], feature[2]));
            });
            cnvnatorPredictMap.put(sample, cnvnatorRegions);
            readKnownCNVRegions(knownCNVFilePath, sample);
        }

        for (Map.Entry<String, Set<Region>> entry : cnvnatorPredictMap.entrySet()) {
            if (!knownCNVMap.containsKey(entry.getKey())) {
                continue;
            }
            Set<Region> natorRegions = entry.getValue();
            Set<Region> knownRegions = knownCNVMap.get(entry.getKey());
            int predictCount = 0;
            int knownCount = knownRegions.stream().mapToInt(region -> {
                if (region.getChr().equals("chr1") || region.getChr().equals("chr4"))
                    return 1;
                else
                    return 0;
            }).sum();
            Map<Region, Set<Region>> beneathRatioMap = new HashMap<>();
            for (Region predictRegion : natorRegions) {
                //                boolean exceedRatio = false;
                for (Region knownRegion : knownRegions) {
                    if (!predictRegion.getChr().equals("chr1")
                        && !predictRegion.getChr().equals("chr4")) {
                        continue; // 只考虑chr1和chr4
                    }
                    if (knownRegion.intersact(predictRegion)) {
                        long predictLength = knownRegion.overlap(predictRegion);
                        if ((double) predictLength
                            / (double) knownRegion.getLength() >= overlapRatio) {
                            predictCount++;
                            //                            exceedRatio = true;
                        } else {
                            if (!beneathRatioMap.containsKey(knownRegion)) {
                                Set<Region> set = new HashSet<>();
                                set.add(predictRegion);
                                beneathRatioMap.put(knownRegion, set);
                            } else {
                                beneathRatioMap.get(knownRegion).add(predictRegion);
                            }
                        }
                    }
                }
                //                if (exceedRatio && beneathRatioMap.containsKey(knownRegion)) {
                //                    beneathRatioMap.remove(knownRegion);
                //                }
            }
            //            predictCount += sumOverlap(beneathRatioMap, overlapRatio);
            if (knownCount != 0) {
                //                System.out.println("CNVnator: " + entry.getKey() + ": " + predictCount
                //                        + ", known CNV regions: " + knownCount);
                System.out.println("CNVnator: " + entry.getKey() + ": " + predictCount
                                   + ", cnvnatorPredictRegions: " + natorRegions.size());
            }
        }
    }

    public static void xhmm(float overlapRatio, String xhmmResultFilePath, String knownCNVFolder) {

        File knownCNVFolderFile = new File(knownCNVFolder);
        String[] knownCNVFileNames = knownCNVFolderFile.list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                if (name.endsWith("tsv"))
                    return true;
                return false;
            }
        });
        List<String> xhmmResultLines = readLines(xhmmResultFilePath);
        xhmmResultLines.remove(0); // drop header
        xhmmResultLines.stream().forEach(line -> {
            String[] feature = line.split("\\s+");
            String[] regionFeature = feature[2].replace(":", "-").split("-");
            Region region = new Region(
                regionFeature[0].length() < 3 ? "chr" + regionFeature[0] : regionFeature[0],
                regionFeature[1], regionFeature[2]);
            String sample = feature[0];
            if (xhmmPredictMap.containsKey(sample)) {
                xhmmPredictMap.get(sample).add(region);
            } else {
                Set<Region> regions = new HashSet<>();
                regions.add(region);
                xhmmPredictMap.put(sample, regions);
            }
        });

        for (String knownCNVFileName : knownCNVFileNames) {
            String[] tmps = knownCNVFileName.replace("_", ".").split("\\."); // split里默认是正则表达式，要按照.分割的话需要转义
            String sample = tmps[tmps.length - 2];
            String knownCNVFilePath = knownCNVFolder + File.separator + knownCNVFileName;
            readKnownCNVRegions(knownCNVFilePath, sample);
        }

        for (Map.Entry<String, Set<Region>> entry : xhmmPredictMap.entrySet()) {
            if (!knownCNVMap.containsKey(entry.getKey())) {
                continue;
            }
            Set<Region> xhmmRegions = entry.getValue();
            Set<Region> knownRegions = knownCNVMap.get(entry.getKey());
            int predictCount = 0;
            int knownCount = knownRegions.stream().mapToInt(region -> {
                if (region.getChr().equals("chr1") || region.getChr().equals("chr4"))
                    return 1;
                else
                    return 0;
            }).sum();
            Map<Region, Set<Region>> beneathRatioMap = new HashMap<>();
            for (Region predictRegion : xhmmRegions) {
                //                boolean exceedRatio = false;
                for (Region knownRegion : knownRegions) {
                    if (!predictRegion.getChr().equals("chr1")
                        && !predictRegion.getChr().equals("chr4")) {
                        continue; // 只考虑chr1和chr4
                    }
                    if (knownRegion.intersact(predictRegion)) {
                        long predictLength = knownRegion.overlap(predictRegion);
                        if ((double) predictLength
                            / (double) knownRegion.getLength() >= overlapRatio) {
                            predictCount++;
                            //                            exceedRatio = true;
                        } else {
                            if (!beneathRatioMap.containsKey(knownRegion)) {
                                Set<Region> set = new HashSet<>();
                                set.add(predictRegion);
                                beneathRatioMap.put(knownRegion, set);
                            } else {
                                beneathRatioMap.get(knownRegion).add(predictRegion);
                            }
                        }
                    }
                }
                //                if (exceedRatio && beneathRatioMap.containsKey(knownRegion)) {
                //                    beneathRatioMap.remove(knownRegion);
                //                }
            }
            //            predictCount += sumOverlap(beneathRatioMap, overlapRatio);
            if (knownCount != 0) {
                //                System.out.println("XHMM: " + entry.getKey() + ": " + predictCount
                //                        + ", known CNV regions: " + knownCount);
                System.out.println("XHMM: " + entry.getKey() + ": " + predictCount
                                   + ", xhmmPredictRegions: " + xhmmRegions.size());
            }
        }

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
    public static int sumOverlap(Map<Region, Set<Region>> beneathRatioMap, float overlapRatio) {
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

    public static void readKnownCNVRegions(String knownCNVFilePath, String sample) {
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
    public static List<String> readLines(String filePath) {
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
