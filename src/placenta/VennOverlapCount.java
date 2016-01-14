package placenta;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import utils.Region;

public class VennOverlapCount {

    private static Set<Region> seqcnvResult   = new HashSet<>();
    private static Set<Region> coniferResult  = new HashSet<>();
    private static Set<Region> cnvnatorResult = new HashSet<>();
    private static Set<Region> cnverResult    = new HashSet<>();
    private static Set<Region> xhmmResult     = new HashSet<>();

    private static Set<Region> knownResult    = new HashSet<>();

    public static void main(String[] args) {
        readAllToolResult();
        List<Set<Region>> list = new ArrayList<>(
            Arrays.asList(seqcnvResult, coniferResult, cnvnatorResult, cnverResult, xhmmResult));
        for (int i = 0; i < list.size(); i++) {
            System.out.println("Count_" + (i + 1) + " = " + list.get(i).size());
        }
        //        System.out.println();
        for (int i = 0; i < list.size() - 1; i++) {
            for (int j = i + 1; j < list.size(); j++) {
                System.out.println("Count_" + (i + 1) + (j + 1) + " = "
                                   + interactAverage(list.get(i), list.get(j)));
            }
        }
        //        System.out.println();
        for (int i = 0; i < (int) Math.ceil(list.size() / 2.0f); i++) {
            for (int j = i + 1; j < list.size() - 1; j++) {
                for (int m = j + 1; m < list.size(); m++)
                    System.out.println("Count_" + (i + 1) + (j + 1) + (m + 1) + " = "
                                       + interactAverage(list.get(i), list.get(j), list.get(m)));
            }
        }
        //        System.out.println();
        for (int i = 0; i < (int) Math.floor(list.size() / 2.0f); i++) {
            for (int j = i + 1; j < (int) Math.ceil(list.size() / 2.0f); j++) {
                for (int m = j + 1; m < list.size() - 1; m++)
                    for (int n = m + 1; n < list.size(); n++)
                        System.out.println("Count_" + (i + 1) + (j + 1) + (m + 1) + (n + 1) + " = "
                                           + interactAverage(list.get(i), list.get(j), list.get(m),
                                               list.get(n)));
            }
        }

        //        System.out.println();
        System.out.println("Count_12345 = " + interactAverage(seqcnvResult, coniferResult,
            cnvnatorResult, cnverResult, xhmmResult));

        System.out.println();
        System.out.println();
        System.out.println();
        readKnownCNVRegions("D:\\placenta_latest\\simulatedRegions.txt");
        System.out.println("SeqCNV vs knownCNV:");
        System.out.println("area1 = " + seqcnvResult.size() + ", area2 = " + knownResult.size()
                           + ", cross.area = " + interactAverage(seqcnvResult, knownResult));

        System.out.println("CoNIFER vs knownCNV:");
        System.out.println("area1 = " + coniferResult.size() + ", area2 = " + knownResult.size()
                           + ", cross.area = " + interactAverage(coniferResult, knownResult));

        System.out.println("CNVnator vs knownCNV:");
        System.out.println("area1 = " + cnvnatorResult.size() + ", area2 = " + knownResult.size()
                           + ", cross.area = " + interactAverage(cnvnatorResult, knownResult));

        System.out.println("CNVer vs knownCNV:");
        System.out.println("area1 = " + cnverResult.size() + ", area2 = " + knownResult.size()
                           + ", cross.area = " + interactAverage(cnverResult, knownResult));

        System.out.println("XHMM vs knownCNV:");
        System.out.println("area1 = " + xhmmResult.size() + ", area2 = " + knownResult.size()
                           + ", cross.area = " + interactAverage(xhmmResult, knownResult));

    }

    /**
     * return average number of regions interacted between set1 and set2, set3, set4, set5.
     * <p>symmetric, e.g. interactAverage(set1, set2, set3, set4, set5) is always equals interactAverage(set3, set4, set2, set1, set5), etc.</p>
     * @param set1
     * @param set2
     * @return
     */
    public static int interactAverage(Set<Region> set1, Set<Region> set2, Set<Region> set3,
                                      Set<Region> set4, Set<Region> set5) {

        List<Set<Region>> list = new ArrayList<>(Arrays.asList(set1, set2, set3, set4, set5));
        list.sort((setItem1, setItem2) -> {
            if (setItem1.size() == setItem2.size())
                return 0;
            if (setItem1.size() < setItem2.size())
                return -1;
            return 1;
        });
        return interact(list.get(0), list.get(1), list.get(2), list.get(3), list.get(4));
        //        int count = 0;
        //        for (int i = 0; i < list.size(); i++) {
        //            for (int j = 0; j < list.size(); j++) {
        //                if (j != i) {
        //                    for (int m = 0; m < list.size(); m++) {
        //                        if (m != i && m != j) {
        //                            for (int n = 0; n < list.size(); n++) {
        //                                if (n != m && n != j && n != i) {
        //                                    for (int x = 0; x < list.size(); x++) {
        //                                        if (x != m && x != n && x != j && x != i) {
        //                                            count += interact(list.get(i), list.get(j), list.get(m),
        //                                                list.get(n), list.get(x));
        //                                        }
        //                                    }
        //                                }
        //                            }
        //                        }
        //                    }
        //                }
        //            }
        //        }
        //        return (int) Math.ceil(count / 120.0f);
    }

    /**
     * return number of regions in set1 which interact with regions in set2, set3, set4, set5.
     * <p><b>Note that the overlap measure is not symmetric, e.g. interact(set1, set2, set3, set4, set5) is not always equals interact(set4, set3, set2, set1, set5) or interact(set1, set4, set3, set2, set5), etc.</b></p>
     * @param set1
     * @param set2
     */
    public static int interact(Set<Region> set1, Set<Region> set2, Set<Region> set3,
                               Set<Region> set4, Set<Region> set5) {
        Set<Region> mergedSet1 = new HashSet<>();
        Set<Region> mergedSet2 = new HashSet<>();
        Set<Region> mergedSet3 = new HashSet<>();

        set1.forEach(region1 -> {
            set2.forEach(region2 -> {
                if (region1.isOverlappedWithType(region2)) {
                    mergedSet1.add(region1);
                }
            });
        });

        if (mergedSet1.isEmpty()) {
            return 0;
        }

        mergedSet1.forEach(region1 -> {
            set3.forEach(region3 -> {
                if (region1.isOverlappedWithType(region3)) {
                    mergedSet2.add(region1);
                }
            });
        });

        if (mergedSet2.isEmpty()) {
            return 0;
        }

        mergedSet2.forEach(region2 -> {
            set4.forEach(region4 -> {
                if (region2.isOverlappedWithType(region4)) {
                    mergedSet3.add(region2);
                }
            });
        });

        if (!mergedSet3.isEmpty()) {
            return mergedSet3.stream().mapToInt(region3 -> {
                for (Region region5 : set5) {
                    if (region3.isOverlappedWithType(region5)) {
                        return 1;
                    }
                }
                return 0;
            }).sum();
        }
        return 0;
    }

    /**
     * return average number of regions interacted between set1 and set2, set3, set4.
     * <p>symmetric, e.g. interactAverage(set1, set2, set3, set4) is always equals interactAverage(set3, set4, set2, set1), etc.</p>
     * @param set1
     * @param set2
     * @return
     */
    public static int interactAverage(Set<Region> set1, Set<Region> set2, Set<Region> set3,
                                      Set<Region> set4) {

        List<Set<Region>> list = new ArrayList<>(Arrays.asList(set1, set2, set3, set4));
        list.sort((setItem1, setItem2) -> {
            if (setItem1.size() == setItem2.size())
                return 0;
            if (setItem1.size() < setItem2.size())
                return -1;
            return 1;
        });
        return interact(list.get(0), list.get(1), list.get(2), list.get(3));
        //        int count = 0;
        //        for (int i = 0; i < list.size(); i++) {
        //            for (int j = 0; j < list.size(); j++) {
        //                if (j != i) {
        //                    for (int m = 0; m < list.size(); m++) {
        //                        if (m != i && m != j) {
        //                            for (int n = 0; n < list.size(); n++) {
        //                                if (n != m && n != j && n != i) {
        //                                    count += interact(list.get(i), list.get(j), list.get(m),
        //                                        list.get(n));
        //                                }
        //                            }
        //                        }
        //                    }
        //                }
        //            }
        //        }
        //        return (int) Math.ceil(count / 24.0f);
    }

    /**
     * return number of regions in set1 which interact with regions in set2, set3, set4.
     * <p><b>Note that the overlap measure is not symmetric, e.g. interact(set1, set2, set3, set4) is not always equals interact(set4, set3, set2, set1) or interact(set1, set4, set3, set2), etc.</b></p>
     * @param set1
     * @param set2
     */
    public static int interact(Set<Region> set1, Set<Region> set2, Set<Region> set3,
                               Set<Region> set4) {
        Set<Region> mergedSet1 = new HashSet<>();
        Set<Region> mergedSet2 = new HashSet<>();

        set1.forEach(region1 -> {
            set2.forEach(region2 -> {
                if (region1.isOverlappedWithType(region2)) {
                    mergedSet1.add(region1);
                }
            });
        });

        if (mergedSet1.isEmpty()) {
            return 0;
        }

        mergedSet1.forEach(region1 -> {
            set3.forEach(region3 -> {
                if (region1.isOverlappedWithType(region3)) {
                    mergedSet2.add(region1);
                }
            });
        });

        if (!mergedSet2.isEmpty()) {
            return mergedSet2.stream().mapToInt(region2 -> {
                for (Region region4 : set4) {
                    if (region2.isOverlappedWithType(region4)) {
                        return 1;
                    }
                }
                return 0;
            }).sum();
        }
        return 0;
    }

    /**
     * return average number of regions interacted between set1 and set2, set3.
     * <p>symmetric, e.g. interactAverage(set1, set2, set3) is always equals interactAverage(set3, set2, set1), etc.</p>
     * @param set1
     * @param set2
     * @return
     */
    public static int interactAverage(Set<Region> set1, Set<Region> set2, Set<Region> set3) {
        List<Set<Region>> list = new ArrayList<>(Arrays.asList(set1, set2, set3));
        list.sort((setItem1, setItem2) -> {
            if (setItem1.size() == setItem2.size())
                return 0;
            if (setItem1.size() < setItem2.size())
                return -1;
            return 1;
        });
        return interact(list.get(0), list.get(1), list.get(2));
        //        return (int) Math.ceil(interact(set1, set2, set3) + interact(set1, set3, set2)
        //                               + interact(set2, set1, set3) + interact(set2, set3, set1)
        //                               + interact(set3, set1, set2) + interact(set3, set2, set1) / 6.0f);
    }

    /**
     * return number of regions in set1 which interact with regions in set2, set3.
     * <p><b>Note that the overlap measure is not symmetric, e.g. interact(set1, set2, set3) is not always equals interact(set3, set2, set1) or interact(set1, set3, set2), etc.</b></p>
     * @param set1
     * @param set2
     */
    public static int interact(Set<Region> set1, Set<Region> set2, Set<Region> set3) {
        Set<Region> mergedSet1 = new HashSet<>();

        set1.forEach(region1 -> {
            set2.forEach(region2 -> {
                if (region1.isOverlappedWithType(region2)) {
                    mergedSet1.add(region1);
                }
            });
        });

        if (!mergedSet1.isEmpty()) {
            return mergedSet1.stream().mapToInt(region1 -> {
                for (Region region2 : set3) {
                    if (region1.isOverlappedWithType(region2)) {
                        return 1;
                    }
                }
                return 0;
            }).sum();
        }
        return 0;
    }

    /**
     * return average number of regions interacted between set1 and set2.
     * <p>symmetric, e.g. interactAverage(set1, set2) is always equals interactAverage(set2, set1).</p>
     * @param set1
     * @param set2
     * @return
     */
    public static int interactAverage(Set<Region> set1, Set<Region> set2) {
        Set<Region> tempSet1 = set1.size() <= set2.size() ? set1 : set2;
        Set<Region> tempSet2 = set1.size() > set2.size() ? set1 : set2;
        return interact(tempSet1, tempSet2);

        //        return (int) Math.ceil(interact(set1, set2) + interact(set2, set1) / 2.0f);
    }

    /**
     * return number of regions in set1 which interact with regions in set2.
     * <p><b>Note that the overlap measure is not symmetric, e.g. interact(set1, set2) is not always equals interact(set2, set1).</b></p>
     * @param set1
     * @param set2
     */
    public static int interact(Set<Region> set1, Set<Region> set2) {
        return set1.stream().mapToInt(region1 -> {
            for (Region region2 : set2) {
                if (region1.isOverlappedWithType(region2)) {
                    return 1;
                }
            }
            return 0;
        }).sum();
    }

    /**
     * read all tool results to memory.
     */
    public static void readAllToolResult() {
        String seqcnvResultPath = "D:\\placenta_resultset\\seqcnv_2stage_2timesPenalty_librarysize\\report\\CNV_report.txt";
        String coniferResultPath = "D:\\placenta_latest\\CoNIFER\\svd_5.txt";
        String cnvnatorResultPath = "D:\\placenta_latest\\CNVnator\\placenta_BAC_predict.txt";
        String cnverResultFolder = "D:\\placenta_latest\\CNVer";
        String xhmmResultPath = "D:\\placenta_latest\\XHMM\\DATA.xcnv";

        List<String> resultLines = readLines(seqcnvResultPath);
        resultLines.forEach(line -> {
            String[] feature = line.split("\\s+");
            feature[0] = feature[0].length() < 3 ? "chr" + feature[0] : feature[0];
            Region predictedRegion = new Region(feature[0], feature[1], feature[2]);
            if (Float.parseFloat(feature[6]) <= 0.6) {
                predictedRegion.setType("LOSS");
            }
            seqcnvResult.add(predictedRegion);
        });

        resultLines = readLines(coniferResultPath);
        resultLines.forEach(line -> {
            String[] feature = line.split("\\s+");
            feature[1] = feature[1].length() < 3 ? "chr" + feature[1] : feature[1];
            Region region = new Region(feature[1], feature[2], feature[3]);
            if (feature[4].trim().equals("del")) {
                region.setType("LOSS");
            }
            coniferResult.add(region);
        });

        resultLines = readLines(cnvnatorResultPath);
        resultLines.forEach(line -> {
            String[] feature = line.split("\\s+");
            String[] cnvFeature = feature[1].replace(":", "-").split("-");
            cnvFeature[0] = cnvFeature[0].length() < 3 ? "chr" + cnvFeature[0] : cnvFeature[0];
            Region predictedRegion = new Region(cnvFeature[0], cnvFeature[1], cnvFeature[2]);
            if (feature[0].trim().equals("deletion")) {
                predictedRegion.setType("LOSS");
            }
            cnvnatorResult.add(predictedRegion);
        });

        File folderFile = new File(cnverResultFolder);
        String[] resultFileNames = folderFile.list((dir, name) -> name.endsWith("cnvs"));
        for (String fileName : resultFileNames) {
            resultLines = readLines(cnverResultFolder + File.separator + fileName);
            resultLines.forEach(line -> {
                String[] feature = line.split("\\s+");
                feature[0] = feature[0].length() < 3 ? "chr" + feature[0] : feature[0];
                Region predictedRegion = new Region(feature[0], feature[1], feature[2]);
                if (feature[3].trim().equals("loss")) {
                    predictedRegion.setType("LOSS");
                }
                cnverResult.add(predictedRegion);
            });
        }

        resultLines = readLines(xhmmResultPath);
        resultLines.remove(0); // drop header
        resultLines.forEach(line -> {
            String[] feature = line.split("\\s+");
            String[] cnvFeature = feature[2].replace(":", "-").split("-");
            cnvFeature[0] = cnvFeature[0].length() < 3 ? "chr" + cnvFeature[0] : cnvFeature[0];
            Region predictRegion = new Region(cnvFeature[0], cnvFeature[1], cnvFeature[2]);
            if (feature[1].trim().equals("DEL")) {
                predictRegion.setType("LOSS");
            }
            xhmmResult.add(predictRegion);
        });

    }

    private static void readKnownCNVRegions(String knownCNVFilePath) {
        if (!knownResult.isEmpty()) {
            return;
        }
        List<String> knownCNVLines = readLines(knownCNVFilePath);
        knownCNVLines.forEach(line -> {
            String[] feature = line.split("\\s+");
            feature[0] = feature[0].length() < 3 ? "chr" + feature[0] : feature[0];
            Region knownRegion = new Region(feature[0], feature[1], feature[2]); // only simulated duplication in placenta dataset.
            knownResult.add(knownRegion);
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
