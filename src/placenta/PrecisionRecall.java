package placenta;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import utils.Region;

/**
 * calculate precision & recall for placenta data. 
 * @author racing
 * @version $Id: PrecisionRecall.java, v 0.1 Dec 2, 2015 19:37:39 PM racing Exp $
 */
public class PrecisionRecall {

    private static final Set<Region> seqcnvPredictRegions   = new HashSet<>();
    private static final Set<Region> coniferPredictRegions  = new HashSet<>();
    private static final Set<Region> cnvnatorPredictRegions = new HashSet<>();
    private static final Set<Region> cnverPredictRegions    = new HashSet<>();
    private static final Set<Region> xhmmPredictRegions     = new HashSet<>();
    private static final Set<Region> knownCNVRegions        = new HashSet<>();

    public static void main(String[] args) {

        System.out.println("Placenta data precision recall analysis:");
        System.out.println("Recall = predicted known CNV events / total known CNV events");
        System.out.println("Precision = unique correctly detected events / total tool CNV events");
        System.out.println();

        float[] overlapRatios = { 0.001f, 0.01f, 0.1f, 0.5f };
        for (float overlapRatio : overlapRatios) {
            // mac
            seqcnv(overlapRatio, "/Users/racing/Desktop/placenta3.5/report/20CNV_report.txt",
                "/Users/racing/Downloads/placenta_latest/simulatedRegions.txt");
            conifer(overlapRatio, "/Users/racing/Downloads/placenta_latest/CoNIFER/svd_5.txt",
                "/Users/racing/Downloads/placenta_latest/simulatedRegions.txt");
            cnvnator(overlapRatio,
                "/Users/racing/Downloads/placenta_latest/CNVnator/placenta_BAC_predict.txt",
                "/Users/racing/Downloads/placenta_latest/simulatedRegions.txt");
            cnver(overlapRatio, "/Users/racing/Downloads/placenta_latest/CNVer",
                "/Users/racing/Downloads/placenta_latest/simulatedRegions.txt");
            xhmm(overlapRatio, "/Users/racing/Downloads/placenta_latest/XHMM/DATA.xcnv",
                "/Users/racing/Downloads/placenta_latest/simulatedRegions.txt");
            // windows
            /*seqcnv(overlapRatio,
            "C:\\Users\\Administrator\\Desktop\\placenta_parameter\\placenta3.5_CNV_report.txt",
            "D:\\placenta_latest\\simulatedRegions.txt");
            conifer(overlapRatio, "D:\\placenta_latest\\CoNIFER\\svd_5.txt",
            "D:\\placenta_latest\\simulatedRegions.txt");
            cnvnator(overlapRatio, "D:\\placenta_latest\\CNVnator\\placenta_BAC_predict.txt",
            "D:\\placenta_latest\\simulatedRegions.txt");
            cnver(overlapRatio, "D:\\placenta_latest\\CNVer",
            "D:\\placenta_latest\\simulatedRegions.txt");
            xhmm(overlapRatio, "D:\\placenta_latest\\XHMM\\DATA.xcnv",
            "D:\\placenta_latest\\simulatedRegions.txt");
            System.out.println();
            System.out.println();
            System.out.println();*/
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
     * @param seqcnvResultPath
     * @param knownCNVPath
     */
    public static void seqcnv(float overlapRatio, String seqcnvResultPath, String knownCNVPath) {
        System.out.println("SeqCNV:");
        System.out.println("Overlap ratio threshold is " + overlapRatio);

        List<String> resultLines = readLines(seqcnvResultPath);
        resultLines.stream().forEach(line -> {
            String[] feature = line.split("\\s+");
            feature[0] = feature[0].length() < 3 ? "chr" + feature[0] : feature[0];
            Region predictedRegion = new Region(feature[0], feature[1], feature[2]);
            if (Float.parseFloat(feature[6]) <= 0.6) {
                predictedRegion.setType("LOSS");
            }
            seqcnvPredictRegions.add(predictedRegion);
        });
        outputPrecisionRecall(overlapRatio, seqcnvPredictRegions, knownCNVPath);
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
     * @param knownCNVPath
     */
    public static void conifer(float overlapRatio, String coniferResultFilePath,
                               String knownCNVPath) {
        System.out.println("CoNIFER:");
        System.out.println("Overlap ratio threshold is " + overlapRatio);

        List<String> resultLines = readLines(coniferResultFilePath);
        resultLines.stream().forEach(line -> {
            String[] feature = line.split("\\s+");
            feature[1] = feature[1].length() < 3 ? "chr" + feature[1] : feature[1];
            Region region = new Region(feature[1], feature[2], feature[3]);
            if (feature[4].trim().equals("del")) {
                region.setType("LOSS");
            }
            coniferPredictRegions.add(region);
        });
        outputPrecisionRecall(overlapRatio, coniferPredictRegions, knownCNVPath);
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
     * @param cnvnatorResultFilePath
     * @param knownCNVPath
     */
    public static void cnvnator(float overlapRatio, String cnvnatorResultFilePath,
                                String knownCNVPath) {
        System.out.println("CNVnator:");
        System.out.println("Overlap ratio threshold is " + overlapRatio);

        List<String> resultLines = readLines(cnvnatorResultFilePath);
        resultLines.stream().forEach(line -> {
            String[] feature = line.split("\\s+");
            String[] cnvFeature = feature[1].replace(":", "-").split("-");
            cnvFeature[0] = cnvFeature[0].length() < 3 ? "chr" + cnvFeature[0] : cnvFeature[0];
            Region predictedRegion = new Region(cnvFeature[0], cnvFeature[1], cnvFeature[2]);
            if (feature[0].trim().equals("deletion")) {
                predictedRegion.setType("LOSS");
            }
            cnvnatorPredictRegions.add(predictedRegion);
        });
        outputPrecisionRecall(overlapRatio, cnvnatorPredictRegions, knownCNVPath);
    }

    /**
     * CNVer的overlap分析，输出precision和recall。
     * <pre>
     * Note:
     * 1. 没有考虑CNVer的CNV region可能会包含多个known CNV region的情况。
     * 2. 如果known CNV region overlap了好几个predict CNV region, 会导致重复计算。
     * 
     * Problem：
     * knownCNV中重复的region要不要去掉？
     * 1. 不去掉，干脆重复计算。
     * 2. 去掉更短的重复region。（待做）
     * </pre>
     * @param overlapRatio
     * @param cnvnatorResultFilePath
     * @param knownCNVPath
     */
    public static void cnver(float overlapRatio, String cnverResultFolder, String knownCNVPath) {
        System.out.println("CNVer:");
        System.out.println("Overlap ratio threshold is " + overlapRatio);

        File folderFile = new File(cnverResultFolder);
        String[] resultFileNames = folderFile.list(new FilenameFilter() {
            @Override
            public boolean accept(File dir, String name) {
                if (name.endsWith("cnvs"))
                    return true;
                return false;
            }
        });
        for (String fileName : resultFileNames) {
            List<String> resultLines = readLines(cnverResultFolder + File.separator + fileName);
            resultLines.stream().forEach(line -> {
                String[] feature = line.split("\\s+");
                feature[0] = feature[0].length() < 3 ? "chr" + feature[0] : feature[0];
                Region predictedRegion = new Region(feature[0], feature[1], feature[2]);
                if (feature[3].trim().equals("loss")) {
                    predictedRegion.setType("LOSS");
                }
                cnverPredictRegions.add(predictedRegion);
            });
        }
        outputPrecisionRecall(overlapRatio, cnverPredictRegions, knownCNVPath);
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
    public static void xhmm(float overlapRatio, String xhmmResultFilePath, String knownCNVPath) {
        System.out.println("XHMM:");
        System.out.println("Overlap ratio threshold is " + overlapRatio);

        List<String> resultLines = readLines(xhmmResultFilePath);
        resultLines.remove(0);
        resultLines.stream().forEach(line -> {
            String[] feature = line.split("\\s+");
            String[] cnvFeature = feature[2].replace(":", "-").split("-");
            cnvFeature[0] = cnvFeature[0].length() < 3 ? "chr" + cnvFeature[0] : cnvFeature[0];
            Region predictRegion = new Region(cnvFeature[0], cnvFeature[1], cnvFeature[2]);
            if (feature[1].trim().equals("DEL")) {
                predictRegion.setType("LOSS");
            }
            xhmmPredictRegions.add(predictRegion);
        });
        outputPrecisionRecall(overlapRatio, xhmmPredictRegions, knownCNVPath);
    }

    /**
     * read known CNV regions and calculate precision & recall.
     * @param overlapRatio
     * @param toolPredictRegions
     * @param knownCNVPath
     */
    public static void outputPrecisionRecall(Float overlapRatio, Set<Region> toolPredictRegions,
                                             String knownCNVPath) {

        Set<Region> uniquePredictedRegions = new HashSet<>();

        int predictedKnownRegions = 0;
        int uniqueCorrectedCNVRegions = 0;
        int knownRegions = 0;
        int toolCNVRegions = 0;

        readKnownCNVRegions(knownCNVPath);

        toolCNVRegions += toolPredictRegions.size();
        knownRegions += knownCNVRegions.size();

        Map<Region, Set<Region>> beneathMap = new HashMap<>(); // overlap没有超过overlapRatio的那些regions，key为knownRegion，value为predictRegions
        Set<Region> predictedKnownRegionSet = new HashSet<>(); // knownRegions中，被predict出来的那些region

        for (Region predictRegion : toolPredictRegions) {
            for (Region knownRegion : knownCNVRegions) {
                if (knownRegion.isOverlappedWithType(predictRegion)) {
                    long predictLength = knownRegion.getOverlapLengthWithType(predictRegion);
                    //                    if ((float) predictLength / (float) knownRegion.getLength() >= overlapRatio) {
                    if ((float) predictLength / (float) predictRegion
                        .getOverlapBaseLengthWithType(knownRegion) >= overlapRatio) {
                        uniquePredictedRegions.add(predictRegion);
                        predictedKnownRegionSet.add(knownRegion);
                    } else {
                        if (beneathMap.containsKey(knownRegion)) {
                            beneathMap.get(knownRegion).add(predictRegion);
                        } else {
                            Set<Region> regions = new HashSet<>();
                            regions.add(predictRegion);
                            beneathMap.put(knownRegion, regions);
                        }
                    }
                }
            }
        }

        for (Map.Entry<Region, Set<Region>> entry : beneathMap.entrySet()) {
            Region knownRegion = entry.getKey();
            Set<Region> predictRegions = entry.getValue();
            long overlapLength = 0;
            List<Region> mergedRegions = new ArrayList<>(); // 计算knownRegion和predictRegions的overlapBaseLength
            mergedRegions.add(knownRegion);
            mergedRegions.addAll(predictRegions);
            Region mergedPredictRegion = Region.mergeOverlappedRegions(mergedRegions);
            for (Region predictRegion : predictRegions) {
                overlapLength += knownRegion.getOverlapLengthWithType(predictRegion);
            }
            if ((double) overlapLength / mergedPredictRegion.getLength() >= overlapRatio) {
                predictedKnownRegionSet.add(knownRegion);
            }
        }
        predictedKnownRegions = predictedKnownRegionSet.size();

        uniqueCorrectedCNVRegions = uniquePredictedRegions.size();

        System.out.println("Correctly detected CNV events: " + predictedKnownRegions
                           + ", Unique correctly detected CNV events: " + uniqueCorrectedCNVRegions
                           + ", Known CNV events: " + knownRegions + ", Tool CNV events: "
                           + toolCNVRegions);
        System.out
            .println(
                "Recall: " + String.format("%.3f", (double) predictedKnownRegions / knownRegions)
                     + ", Precision: "
                     + String.format("%.3f", (double) uniqueCorrectedCNVRegions / toolCNVRegions));
        System.out.println();
    }

    /**
     * read known CNV regions and calculate precision & recall.
     * @param overlapRatio
     * @param toolPredictRegions
     * @param knownCNVPath
     */
    public static void outputPrecisionRecallBreakDancer(Float overlapRatio,
                                                        Set<Region> toolPredictRegions,
                                                        String knownCNVPath) {
        overlapRatio = 0.5f;
        Set<Region> uniquePredictedRegions = new HashSet<>();

        int predictedKnownRegions = 0;
        int uniqueCorrectedCNVRegions = 0;
        int knownRegions = 0;
        int toolCNVRegions = 0;

        readKnownCNVRegions(knownCNVPath);

        toolCNVRegions += toolPredictRegions.size();
        knownRegions += knownCNVRegions.size();

        Map<Region, Set<Region>> beneathMap = new HashMap<>(); // overlap没有超过overlapRatio的那些regions，key为knownRegion，value为predictRegions
        Set<Region> predictedKnownRegionSet = new HashSet<>(); // knownRegions中，被predict出来的那些region

        for (Region predictRegion : toolPredictRegions) {
            for (Region knownRegion : knownCNVRegions) {
                if (knownRegion.isOverlappedWithType(predictRegion)) {
                    long predictLength = knownRegion.getOverlapLengthWithType(predictRegion);

                    if (predictRegion.getType().toString().equals("GAIN")) {
                        if (predictLength > 0) {
                            uniquePredictedRegions.add(predictRegion);
                            predictedKnownRegionSet.add(knownRegion);
                        }
                    } else {
                        if ((float) predictLength
                            / (float) predictRegion.getLength() >= overlapRatio) {
                            if ((float) predictLength
                                / (float) knownRegion.getLength() >= overlapRatio) {
                                uniquePredictedRegions.add(predictRegion);
                                predictedKnownRegionSet.add(knownRegion);
                            } else {
                                if (beneathMap.containsKey(knownRegion)) {
                                    beneathMap.get(knownRegion).add(predictRegion);
                                } else {
                                    Set<Region> regions = new HashSet<>();
                                    regions.add(predictRegion);
                                    beneathMap.put(knownRegion, regions);
                                }
                            }
                        }
                    }
                }
            }
        }

        for (Map.Entry<Region, Set<Region>> entry : beneathMap.entrySet()) {
            Region knownRegion = entry.getKey();
            Set<Region> predictRegions = entry.getValue();
            long overlapLength = 0;
            List<Region> mergedRegions = new ArrayList<>(); // 计算knownRegion和predictRegions的overlapBaseLength
            mergedRegions.add(knownRegion);
            mergedRegions.addAll(predictRegions);
            Region mergedPredictRegion = Region.mergeOverlappedRegions(mergedRegions);
            for (Region predictRegion : predictRegions) {
                overlapLength += knownRegion.getOverlapLengthWithType(predictRegion);
            }
            if ((double) overlapLength / mergedPredictRegion.getLength() >= overlapRatio) {
                predictedKnownRegionSet.add(knownRegion);
            }
        }
        predictedKnownRegions = predictedKnownRegionSet.size();

        uniqueCorrectedCNVRegions = uniquePredictedRegions.size();

        System.out.println("Correctly detected CNV events: " + predictedKnownRegions
                           + ", Unique correctly detected CNV events: " + uniqueCorrectedCNVRegions
                           + ", Known CNV events: " + knownRegions + ", Tool CNV events: "
                           + toolCNVRegions);
        System.out
            .println(
                "Recall: " + String.format("%.3f", (double) predictedKnownRegions / knownRegions)
                     + ", Precision: "
                     + String.format("%.3f", (double) uniqueCorrectedCNVRegions / toolCNVRegions));
        System.out.println();
    }

    private static void readKnownCNVRegions(String knownCNVFilePath) {
        if (!knownCNVRegions.isEmpty()) {
            return;
        }
        List<String> knownCNVLines = readLines(knownCNVFilePath);
        knownCNVLines.stream().forEach(line -> {
            String[] feature = line.split("\\s+");
            feature[0] = feature[0].length() < 3 ? "chr" + feature[0] : feature[0];
            Region knownRegion = new Region(feature[0], feature[1], feature[2]); // only simulated duplication in placenta dataset.
            knownCNVRegions.add(knownRegion);
        });
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
