package wes;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
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

        float[] overlapRatios = { 0.01f, 0.1f, 0.5f };
        for (float overlapRatio : overlapRatios) {
            //            seqcnv(overlapRatio, "C:\\Users\\Administrator\\Desktop\\chr1_chr4_result\\SeqCNV",
            //                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv"); // chr1_chr4
            //            seqcnv(overlapRatio, "E:\\SeqCNV_previous_mergedBed",
            //                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv"); // mergedBed
            //            seqcnv(overlapRatio, "E:\\SeqCNV_previous_nonMergedBed",
            //                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv"); // previous_librarysize
            //
            //            seqcnv(overlapRatio, "D:\\seqcnv_2stage_weschr1chr4",
            //                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv"); // 2stage
            //            seqcnv(overlapRatio, "D:\\seqcnv_2stage_weschr1chr4_merged",
            //                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv"); // 2stage_merged
            //            seqcnv(overlapRatio, "D:\\seqcnv_2stage_2timesPenalty_weschr1chr4",
            //                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv"); // 2stage_2timesPenalty
            seqcnv(overlapRatio, "D:\\seqcnv_2stage_2timesPenalty_weschr1chr4_merged",
                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv"); // 2stage_2timesPenalty_merged
            //            seqcnv(overlapRatio, "C:\\Users\\Administrator\\Desktop\\seqcnv_chr1",
            //                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv"); // only chr1
            //            seqcnv(overlapRatio, "C:\\Users\\Administrator\\Desktop\\seqcnv_penalty_chr4",
            //                "C:\\Users\\Administrator\\Desktop\\known_cnv\\dgv_cnv"); // only chr4
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
        //        System.out.println("Overlap ratio threshold is " + overlapRatio);

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
                if (analysedChromosomes.contains(feature[0])) {
                    Region predictedRegion = new Region(feature[0], feature[1], feature[2]);
                    if (Float.parseFloat(feature[6]) <= 0.6) {
                        predictedRegion.setType("LOSS");
                    }
                    seqcnvRegions.add(predictedRegion);
                }
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
        //        System.out.println("Overlap ratio threshold is " + overlapRatio);

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
            if (feature[4].trim().equals("del")) {
                region.setType("LOSS");
            }
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
        //        System.out.println("Overlap ratio threshold is " + overlapRatio);

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
                String[] lineFeature = line.split("\\s+");
                String[] feature = lineFeature[1].replace(":", "-").split("-");
                feature[0] = feature[0].length() < 3 ? "chr" + feature[0] : feature[0];
                if (analysedChromosomes.contains(feature[0])) {
                    Region predictedRegion = new Region(feature[0], feature[1], feature[2]);
                    if (lineFeature[0].trim().equals("deletion")) {
                        predictedRegion.setType("LOSS");
                    }
                    cnvnatorRegions.add(predictedRegion);
                }
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
        //        System.out.println("Overlap ratio threshold is " + overlapRatio);

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
            if (feature[1].trim().equals("DEL")) {
                region.setType("LOSS");
            }
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
        //        System.out.println("Overlap ratio threshold is " + overlapRatio);

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
                if (analysedChromosomes.contains(feature[0])) {
                    Region predictedRegion = new Region(feature[0], feature[1], feature[2]);
                    if (Integer.parseInt(feature[6].trim()) < 0) {
                        predictedRegion.setType("LOSS");
                    }
                    excavatorRegions.add(predictedRegion);
                }
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
        outputPrecisionRecallBreakDancer(overlapRatio, toolPredictMap, knownCNVFolder);

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
                    if (knownRegion.isOverlappedWithType(predictRegion)) {
                        long predictLength = knownRegion.getOverlapLengthWithType(predictRegion);
                        if ((double) predictLength / (double) predictRegion
                            .getOverlapBaseLengthWithType(knownRegion) >= overlapRatio) {
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
            }

            int extraPredictedKnown = 0;
            Set<Region> extraCorrectlyPredictedRegions = new HashSet<>();
            Map<Region, Set<Region>> sampleBeneathMap = beneathMap.get(sample);
            if (sampleBeneathMap != null && !sampleBeneathMap.isEmpty()) {
                for (Map.Entry<Region, Set<Region>> beneathEntry : sampleBeneathMap.entrySet()) {
                    Region knownRegion = beneathEntry.getKey();
                    Set<Region> predictRegions = beneathEntry.getValue();
                    List<Region> mergedRegions = new ArrayList<>(); // 计算knownRegion和predictRegions的overlapBaseLength
                    mergedRegions.add(knownRegion);
                    mergedRegions.addAll(predictRegions);
                    Region mergedPredictRegion = Region.mergeOverlappedRegions(mergedRegions);

                    long overlapLength = 0;
                    for (Region predictRegion : predictRegions) {
                        overlapLength += knownRegion.getOverlapLengthWithType(predictRegion);
                    }

                    if ((double) overlapLength / mergedPredictRegion.getLength() >= overlapRatio) {
                        extraPredictedKnown++;
                        extraCorrectlyPredictedRegions.addAll(predictRegions);
                    }
                }
            }

            float sampleRecall = (float) (samplePredictedKnownRegions.size() + extraPredictedKnown)
                                 / sampleKnownCNVRegions.size();
            float samplePrecision = (float) (sampleCorrectedCNVRegions.size()
                                             + extraCorrectlyPredictedRegions.size())
                                    / sampleToolCNVRegions.size();
            //            System.out.println(sample + ":" + String.format("%.3f", sampleRecall) + ", "
            //                               + String.format("%.3f", samplePrecision));
            preRecMap.put(sample, new Pair<>(sampleRecall, samplePrecision));

        }

        sampleSet.addAll(excludedSamples);
        sampleSet.removeAll(excludedSamples);

        int analysedSampleCount = sampleSet.size();
        double avgRecall = preRecMap.entrySet().stream().mapToDouble(entry -> {
            return entry.getValue().getFirst();
        }).sum() / analysedSampleCount;
        double avgPrecision = preRecMap.entrySet().stream().mapToDouble(entry -> {
            return entry.getValue().getSecond();
        }).sum() / analysedSampleCount;

        //        BigDecimal avgRecallRounded = new BigDecimal((double) avgRecall);
        //        avgRecallRounded = avgRecallRounded.setScale(2, RoundingMode.HALF_UP);
        //        avgRecall = avgRecallRounded.floatValue();
        //
        //        BigDecimal avgPrecisionRounded = new BigDecimal((double) avgPrecision);
        //        avgPrecisionRounded = avgPrecisionRounded.setScale(2, RoundingMode.HALF_UP);
        //        avgPrecision = avgPrecisionRounded.floatValue();

        System.out.println(String.format("%16s", "OverlapRatio " + overlapRatio) + ": Recall: "
                           + String.format("%.2f", avgRecall) + ", Precision: "
                           + String.format("%.2f", avgPrecision) + " (average)");

    }

    /**
     * read known CNV regions and calculate precision & recall using BreakDancer standard:
     * <p>
     * A loss would not be considered as detected unless its overlapped region with a 
     * predicted loss exceeds 50%, mutually. Since it’s more difficult to detect, 
     * a gain was considered as detected once it overlapped a predicted gain.
     * </p>
     * If the overlap length exceeds 1bp(for copy GAIN), the region is correctly detected.
     * If the overlap length doesn't exceed 50% of the predict Region(for copy LOSS), discard it.
     * If the overlap length exceed 50% of both the predict Region and the known Region, then the region is correctly detected.
     * @param overlapRatio
     * @param toolPredictMap
     * @param knownCNVFolder
     */
    public static void outputPrecisionRecallBreakDancer(Float overlapRatio,
                                                        Map<String, Set<Region>> toolPredictMap,
                                                        String knownCNVFolder) {
        overlapRatio = 0.5f;
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
                    if (knownRegion.isOverlappedWithType(predictRegion)) {

                        long predictLength = knownRegion.getOverlapLengthWithType(predictRegion);

                        if (predictRegion.getType().toString().equals("GAIN")) { // the CNV region is copy gain region
                            if (predictLength > 0) {
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
                        } else { // the CNV region is copy loss region
                            if ((double) predictLength
                                / (double) predictRegion.getLength() >= overlapRatio) {
                                if ((double) predictLength
                                    / (double) knownRegion.getLength() >= overlapRatio) {
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
                                } else {
                                    if (beneathMap.containsKey(sample)) {
                                        Map<Region, Set<Region>> beneathMapItem = beneathMap
                                            .get(sample);
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
                    }
                }
            }

            int extraPredictedKnown = 0;
            Set<Region> extraCorrectlyPredictedRegions = new HashSet<>();
            Map<Region, Set<Region>> sampleBeneathMap = beneathMap.get(sample);
            if (sampleBeneathMap != null && !sampleBeneathMap.isEmpty()) {
                for (Map.Entry<Region, Set<Region>> beneathEntry : sampleBeneathMap.entrySet()) {
                    Region knownRegion = beneathEntry.getKey();
                    Set<Region> predictRegions = beneathEntry.getValue();
                    List<Region> mergedRegions = new ArrayList<>(); // 计算knownRegion和predictRegions的overlapBaseLength
                    mergedRegions.add(knownRegion);
                    mergedRegions.addAll(predictRegions);
                    Region mergedPredictRegion = Region.mergeOverlappedRegions(mergedRegions);

                    long overlapLength = 0;
                    for (Region predictRegion : predictRegions) {
                        overlapLength += knownRegion.getOverlapLengthWithType(predictRegion);
                    }
                    if ((double) overlapLength / mergedPredictRegion.getLength() >= overlapRatio) {
                        extraPredictedKnown++;
                        extraCorrectlyPredictedRegions.addAll(predictRegions);
                    }
                }
            }

            float sampleRecall = (float) (samplePredictedKnownRegions.size() + extraPredictedKnown)
                                 / sampleKnownCNVRegions.size();
            float samplePrecision = (float) (sampleCorrectedCNVRegions.size()
                                             + extraCorrectlyPredictedRegions.size())
                                    / sampleToolCNVRegions.size();
            //            System.out.println(sample + ":" + String.format("%.2f", sampleRecall) + ", "
            //                               + String.format("%.2f", samplePrecision));
            preRecMap.put(sample, new Pair<>(sampleRecall, samplePrecision));

        }

        sampleSet.addAll(excludedSamples);
        sampleSet.removeAll(excludedSamples);

        int analysedSampleCount = sampleSet.size();
        double avgRecall = preRecMap.entrySet().stream().mapToDouble(entry -> {
            return entry.getValue().getFirst();
        }).sum() / analysedSampleCount;
        double avgPrecision = preRecMap.entrySet().stream().mapToDouble(entry -> {
            return entry.getValue().getSecond();
        }).sum() / analysedSampleCount;

        BigDecimal avgRecallRounded = new BigDecimal((double) avgRecall);
        avgRecallRounded = avgRecallRounded.setScale(2, RoundingMode.HALF_UP);
        avgRecall = avgRecallRounded.floatValue();

        BigDecimal avgPrecisionRounded = new BigDecimal((double) avgPrecision);
        avgPrecisionRounded = avgPrecisionRounded.setScale(2, RoundingMode.HALF_UP);
        avgPrecision = avgPrecisionRounded.floatValue();

        System.out.println(
            String.format("%16s", "BreakDancer") + ": Recall: " + String.format("%.2f", avgRecall)
                           + ", Precision: " + String.format("%.2f", avgPrecision) + " (average)");

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
            if (analysedChromosomes.contains(feature[0])) {
                Region knownRegion = new Region(feature[0], feature[1], feature[2]);
                if (feature[10].trim().equals("loss")) {
                    knownRegion.setType("LOSS");
                }
                knownCNVRegions.add(knownRegion);
            }
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
