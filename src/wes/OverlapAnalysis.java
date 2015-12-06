package wes;

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

public class OverlapAnalysis {

    private static final Map<String, Set<Region>>  knownCNVMap        = new HashMap<>();
    private static final Map<String, List<Region>> seqcnvPredictMap   = new HashMap<>();
    private static final Map<String, List<Region>> coniferPredictMap  = new HashMap<>();
    private static final Map<String, List<Region>> cnvnatorPredictMap = new HashMap<>();
    private static final Map<String, List<Region>> xhmmPredictMap     = new HashMap<>();

    public static void main(String[] args) {
        float overlapRatio = 0.1f;
        seqcnv(overlapRatio, "/Users/racing/Downloads/seqcnv_programs/wes",
            "/Users/racing/Desktop/dgv_cnv");
        conifer(overlapRatio, "/Users/racing/Lab/CNV_WES/CoNIFER/output/WES_CoNIFER_Result.txt",
            "/Users/racing/Desktop/dgv_cnv");
        cnvnator(overlapRatio);
        xhmm(overlapRatio);
    }

    /**
     * SeqCNV的overlap分析。考虑了SeqCNV的CNV region可能会包含多个known CNV region的情况。
     * @param overlapRatio
     * @param seqcnvResultFolder
     * @param knownCNVFolder
     */
    public static void seqcnv(float overlapRatio, String seqcnvResultFolder,
                              String knownCNVFolder) {

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
            String knownCNVFilePath = knownCNVFolder + File.separator + "dbvar_known_cnv_" + sample
                                      + ".tsv";
            List<String> seqcnvResultLines = readLines(resultFilePath);
            List<Region> seqcnvRegions = new ArrayList<>();
            seqcnvResultLines.stream().forEach(line -> {
                String[] feature = line.split("\\s+");
                seqcnvRegions
                    .add(new Region(feature[0].length() < 3 ? "chr" + feature[0] : feature[0],
                        feature[1], feature[2]));
            });
            seqcnvPredictMap.put(sample, seqcnvRegions);
            readKnownCNVRegions(knownCNVFilePath, sample);
        }
        for (Map.Entry<String, List<Region>> entry : seqcnvPredictMap.entrySet()) {
            List<Region> seqRegions = entry.getValue();
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
            for (Region knownRegion : knownRegions) {
                boolean exceedRatio = false;
                for (Region predictRegion : seqRegions) {
                    if (knownRegion.intersact(predictRegion)) {
                        long predictLength = knownRegion.overlap(predictRegion);
                        if ((float) predictLength
                            / (float) knownRegion.getLength() >= overlapRatio) {
                            predictCount++;
                            exceedRatio = true;
                            //                                System.out.println(knownRegion + " " + predictRegion);
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
                if (exceedRatio && beneathRatioMap.containsKey(knownRegion)) {
                    beneathRatioMap.remove(knownRegion);
                }
            }
            predictCount += sumOverlap(beneathRatioMap, overlapRatio);
            System.out.println("SeqCNV: " + entry.getKey() + ": " + predictCount
                               + ", known CNV regions: " + knownCount);
        }
    }

    /**
     * CoNIFER的overlap分析。没有考虑CoNIFER的CNV region可能会包含多个known CNV region的情况。
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
        for (String knownCNVFileName : knownCNVFileNames) {
            String[] tmps = knownCNVFileName.replace("_", ".").split("\\."); // split里默认是正则表达式，要按照.分割的话需要转义
            String sample = tmps[tmps.length - 2];
            String knownCNVFilePath = knownCNVFolder + File.separator + knownCNVFileName;
            List<String> coniferResultLines = readLines(coniferResultFilePath);
            coniferResultLines.stream().forEach(line -> {
                String[] feature = line.split("\\s+");
                Region region = new Region(
                    feature[1].length() < 3 ? "chr" + feature[1] : feature[1], feature[2],
                    feature[3]);
                String coniferSample = feature[0];
                if (coniferPredictMap.containsKey(coniferSample)) {
                    coniferPredictMap.get(coniferSample).add(region);
                } else {
                    List<Region> regions = new ArrayList<>();
                    regions.add(region);
                    coniferPredictMap.put(coniferSample, regions);
                }
            });
            readKnownCNVRegions(knownCNVFilePath, sample);
        }

        for (Map.Entry<String, List<Region>> entry : coniferPredictMap.entrySet()) {
            if (!knownCNVMap.containsKey(entry.getKey())) {
                continue;
            }
            System.out.println(entry.getKey());
            List<Region> conRegions = entry.getValue();
            Set<Region> knownRegions = knownCNVMap.get(entry.getKey());
            int predictCount = 0;
            int knownCount = knownRegions.stream().mapToInt(region -> {
                if (region.getChr().equals("chr1") || region.getChr().equals("chr4"))
                    return 1;
                else
                    return 0;
            }).sum();
            Map<Region, Set<Region>> beneathRatioMap = new HashMap<>();
            for (Region knownRegion : knownRegions) {
                boolean exceedRatio = false;
                for (Region predictRegion : conRegions) {
                    if (!predictRegion.getChr().equals("chr1")
                        && !predictRegion.getChr().equals("chr4")) {
                        continue; // 只考虑chr1和chr4
                    }
                    if (knownRegion.intersact(predictRegion)) {
                        long predictLength = knownRegion.overlap(predictRegion);
                        if ((double) predictLength
                            / (double) knownRegion.getLength() >= overlapRatio) {
                            predictCount++;
                            exceedRatio = true;
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
                if (exceedRatio && beneathRatioMap.containsKey(knownRegion)) {
                    beneathRatioMap.remove(knownRegion);
                }
            }
            predictCount += sumOverlap(beneathRatioMap, overlapRatio);
            if (knownCount != 0) {
                System.out.println("CoNIFER: " + entry.getKey() + ": " + predictCount
                                   + ", known CNV regions: " + knownCount);
            }
        }
    }

    public static void cnvnator(float overlapRatio) {
        String resultFilePath = "/Users/racing/Lab/CNV_WES/CNVnator/output/NA10847/NA10847.call";
        String knownCNVFilePath = "/Users/racing/Desktop/dgv_cnv/dbvar_known_cnv_NA10847.tsv";
        List<String> cnvnatorResultLines = readLines(resultFilePath);
        List<Region> cnvnatorRegions = new ArrayList<>();
        cnvnatorResultLines.stream().forEach(line -> {
            String[] feature = line.split("\\s+")[1].replace(":", "-").split("-");
            cnvnatorRegions.add(new Region(
                feature[0].length() < 3 ? "chr" + feature[0] : feature[0], feature[1], feature[2]));
        });
        cnvnatorPredictMap.put("NA10847", cnvnatorRegions);
        readKnownCNVRegions(knownCNVFilePath, "NA10847");

        for (Map.Entry<String, List<Region>> entry : cnvnatorPredictMap.entrySet()) {
            if (!knownCNVMap.containsKey(entry.getKey())) {
                continue;
            }
            List<Region> natorRegions = entry.getValue();
            Set<Region> knownRegions = knownCNVMap.get(entry.getKey());
            int predictCount = 0;
            int knownCount = knownRegions.stream().mapToInt(region -> {
                if (region.getChr().equals("chr1") || region.getChr().equals("chr4"))
                    return 1;
                else
                    return 0;
            }).sum();
            Map<Region, Set<Region>> beneathRatioMap = new HashMap<>();
            for (Region knownRegion : knownRegions) {
                boolean exceedRatio = false;
                for (Region predictRegion : natorRegions) {
                    if (!predictRegion.getChr().equals("chr1")
                        && !predictRegion.getChr().equals("chr4")) {
                        continue; // 只考虑chr1和chr4
                    }
                    if (knownRegion.intersact(predictRegion)) {
                        long predictLength = knownRegion.overlap(predictRegion);
                        if ((double) predictLength
                            / (double) knownRegion.getLength() >= overlapRatio) {
                            predictCount++;
                            exceedRatio = true;
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
                if (exceedRatio && beneathRatioMap.containsKey(knownRegion)) {
                    beneathRatioMap.remove(knownRegion);
                }
            }
            predictCount += sumOverlap(beneathRatioMap, overlapRatio);
            if (knownCount != 0) {
                System.out.println("CNVnator: " + entry.getKey() + ": " + predictCount
                                   + ", known CNV regions: " + knownCount);
            }
        }
    }

    public static void xhmm(float overlapRatio) {
        String resultFilePath = "/Users/racing/Lab/CNV_WES/XHMM/DATA.xcnv";
        String knownCNVFilePath = "/Users/racing/Desktop/dgv_cnv/dbvar_known_cnv_NA10847.tsv";
        List<String> xhmmResultLines = readLines(resultFilePath);
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
                List<Region> regions = new ArrayList<>();
                regions.add(region);
                xhmmPredictMap.put(sample, regions);
            }
        });
        readKnownCNVRegions(knownCNVFilePath, "NA10847");
        for (Map.Entry<String, List<Region>> entry : xhmmPredictMap.entrySet()) {
            if (!knownCNVMap.containsKey(entry.getKey())) {
                continue;
            }
            List<Region> xhmmRegions = entry.getValue();
            Set<Region> knownRegions = knownCNVMap.get(entry.getKey());
            int predictCount = 0;
            int knownCount = knownRegions.stream().mapToInt(region -> {
                if (region.getChr().equals("chr1") || region.getChr().equals("chr4"))
                    return 1;
                else
                    return 0;
            }).sum();
            Map<Region, Set<Region>> beneathRatioMap = new HashMap<>();
            for (Region knownRegion : knownRegions) {
                boolean exceedRatio = false;
                for (Region predictRegion : xhmmRegions) {
                    if (!predictRegion.getChr().equals("chr1")
                        && !predictRegion.getChr().equals("chr4")) {
                        continue; // 只考虑chr1和chr4
                    }
                    if (knownRegion.intersact(predictRegion)) {
                        long predictLength = knownRegion.overlap(predictRegion);
                        if ((double) predictLength
                            / (double) knownRegion.getLength() >= overlapRatio) {
                            predictCount++;
                            exceedRatio = true;
                            System.out.println(knownRegion + " " + predictRegion);
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
                if (exceedRatio && beneathRatioMap.containsKey(knownRegion)) {
                    beneathRatioMap.remove(knownRegion);
                }
            }
            predictCount += sumOverlap(beneathRatioMap, overlapRatio);
            if (knownCount != 0) {
                System.out.println("XHMM: " + entry.getKey() + ": " + predictCount
                                   + ", known CNV regions: " + knownCount);
            }
        }

    }

    /**
     * sum the overlap between the knownCNVRegion and predictCNVRegion, return the num of knownCNVRegions which was overlapped by predictCNVRegions exceed the overlapRatio.
     * @param beneathRatioMap key is knownCNVRegion, value is a set of predictCNVRegions which overlap with the key but beneath the overlapRatio. 
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
        Set<Region> knownCNVRegions = new HashSet<>();
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
