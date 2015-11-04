package cnvstatistic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import utils.Region;

/**
 * 统计各个工具识别了多少cnv region，overlap了多少cnv region等。
 * @author racing
 * @version $Id: Statistic.java, v 0.1 Nov 2, 2015 9:37:39 PM racing Exp $
 */
public class Statistic {
    private static List<Region>         simulatedRegionList = new ArrayList<Region>();
    private static List<Region>         cnverRegionList     = new ArrayList<Region>();
    private static List<Region>         seqcnvRegionList    = new ArrayList<Region>();
    private static List<Region>         cnvnatorRegionList  = new ArrayList<Region>();
    private static List<Region>         coniferRegionList   = new ArrayList<Region>();
    private static Map<Region, Integer> cnverRegionMap      = new HashMap<Region, Integer>();
    private static Map<Region, Integer> seqcnvRegionMap     = new HashMap<Region, Integer>();
    private static Map<Region, Integer> cnvnatorRegionMap   = new HashMap<Region, Integer>();
    private static Map<Region, Integer> coniferRegionMap    = new HashMap<Region, Integer>();

    public static void main(String[] args) {
        readSimulatedRegions("/Users/racing/Lab/(2015-09-22)CNV/simulatedRegions.txt");
        //        System.out.println(list);
        statisticCNVer("/Users/racing/Lab/(2015-09-22)CNV/CNVer");
        statisticCNVnator("/Users/racing/Lab/(2015-09-22)CNV/CNVnator/placenta_BAC_predict.txt");
        statisticConifer("/Users/racing/Lab/(2015-09-22)CNV/Conifer/calls.txt");
        statisticSeqcnv("/Users/racing/Lab/(2015-09-22)CNV/SeqCNV/CNV_report.txt");

        System.out
            .println("CNVer & CNVnator : " + intersect(cnverRegionList, cnvnatorRegionList).size());
        System.out
            .println("CNVer & CoNIFER : " + intersect(cnverRegionList, coniferRegionList).size());
        System.out
            .println("CNVer & SeqCNV : " + intersect(cnverRegionList, seqcnvRegionList).size());
        System.out.println(
            "CoNIFER & CNVnator : " + intersect(coniferRegionList, cnvnatorRegionList).size());
        System.out
            .println("CoNIFER & SeqCNV : " + intersect(coniferRegionList, seqcnvRegionList).size());
        System.out.println(
            "CNVnator & SeqCNV : " + intersect(cnvnatorRegionList, seqcnvRegionList).size());
    }

    public static boolean judgeInSimulatedRegions(Region region) {
        for (Region item : simulatedRegionList) {
            if (region.getChr().equals(item.getChr())) {
                if (region.getStart() >= item.getStart() && region.getStart() <= item.getEnd()
                    || region.getEnd() >= item.getStart() && region.getEnd() <= item.getEnd()) {
                    return true;
                }
                if (region.getStart() <= item.getStart() && region.getEnd() >= item.getEnd()) {
                    return true;
                }
            }
        }
        return false;
    }

    public static void readSimulatedRegions(String path) {
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String str = null;
            while ((str = br.readLine()) != null) {
                String[] features = str.split("\\s+");
                simulatedRegionList.add(new Region(features[0], Integer.parseInt(features[1]),
                    Integer.parseInt(features[2])));
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void statisticCNVer(String dir) {
        int sum = 0;
        int regionsInSimulated = 0;

        File[] files = new File(dir).listFiles();
        for (File file : files) {
            if (file.getName().endsWith("cnvs")) {
                try (BufferedReader br = new BufferedReader(new FileReader(file))) {
                    String str = null;
                    while ((str = br.readLine()) != null) {
                        String[] features = str.split("\\s+");
                        sum++;
                        Region region = new Region(features[0], Integer.parseInt(features[1]),
                            Integer.parseInt(features[2]));
                        if (judgeInSimulatedRegions(region)) {
                            cnverRegionList.add(region);
                            regionsInSimulated++;
                        }
                    }
                } catch (FileNotFoundException e) {
                    e.printStackTrace();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
        System.out.println("CNVer count:" + sum + ", regionsInSimulated:" + regionsInSimulated);
        Collections.sort(cnverRegionList);
        for (Region region : simulatedRegionList) {
            cnverRegionMap.put(region, 0);
        }
        for (Region region : cnverRegionList) {
            Region ref = null;
            for (Map.Entry<Region, Integer> entry : cnverRegionMap.entrySet()) {
                if (entry.getKey().intersact(region)) {
                    ref = entry.getKey();
                    break;
                }
            }
            if (ref != null) {
                cnverRegionMap.put(ref, cnverRegionMap.get(ref) + 1);
            }
        }
        for (Map.Entry<Region, Integer> entry : cnverRegionMap.entrySet()) {
            System.out.println(entry.getKey() + "   " + entry.getValue());
        }
        System.out.println();
        System.out.println();
        System.out.println();

        //        System.out.println(cnverRegionList);
    }

    public static void statisticCNVnator(String path) {
        int sum = 0;
        int regionsInSimulated = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String str = null;
            while ((str = br.readLine()) != null) {
                String[] features = str.split("\\s+")[1].replace(":", " ").replace("-", " ")
                    .split("\\s+");
                sum++;
                Region region = new Region("chr" + features[0], Integer.parseInt(features[1]),
                    Integer.parseInt(features[2]));
                if (judgeInSimulatedRegions(region)) {
                    cnvnatorRegionList.add(region);
                    regionsInSimulated++;
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("CNVnator count:" + sum + ", regionsInSimulated:" + regionsInSimulated);
        Collections.sort(cnvnatorRegionList);
        for (Region region : simulatedRegionList) {
            cnvnatorRegionMap.put(region, 0);
        }
        for (Region region : cnvnatorRegionList) {
            Region ref = null;
            for (Map.Entry<Region, Integer> entry : cnvnatorRegionMap.entrySet()) {
                if (entry.getKey().intersact(region)) {
                    ref = entry.getKey();
                    break;
                }
            }
            if (ref != null) {
                cnvnatorRegionMap.put(ref, cnvnatorRegionMap.get(ref) + 1);
            }
        }
        for (Map.Entry<Region, Integer> entry : cnvnatorRegionMap.entrySet()) {
            System.out.println(entry.getKey() + "   " + entry.getValue());
        }
        System.out.println();
        System.out.println();
        System.out.println();
        //        System.out.println(cnvnatorRegionList);
    }

    public static void statisticConifer(String path) {
        int sum = 0;
        int regionsInSimulated = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String str = null;
            br.readLine();
            while ((str = br.readLine()) != null) {
                String[] features = str.split("\\s+");
                sum++;
                Region region = new Region(features[1], Integer.parseInt(features[2]),
                    Integer.parseInt(features[3]));
                if (judgeInSimulatedRegions(region)) {
                    coniferRegionList.add(region);
                    regionsInSimulated++;
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("CoNIFER count:" + sum + ", regionsInSimulated:" + regionsInSimulated);
        Collections.sort(coniferRegionList);

        for (Region region : simulatedRegionList) {
            coniferRegionMap.put(region, 0);
        }
        for (Region region : coniferRegionList) {
            Region ref = null;
            for (Map.Entry<Region, Integer> entry : coniferRegionMap.entrySet()) {
                if (entry.getKey().intersact(region)) {
                    ref = entry.getKey();
                    break;
                }
            }
            if (ref != null) {
                coniferRegionMap.put(ref, coniferRegionMap.get(ref) + 1);
            }
        }
        for (Map.Entry<Region, Integer> entry : coniferRegionMap.entrySet()) {
            System.out.println(entry.getKey() + "   " + entry.getValue());
        }
        System.out.println();
        System.out.println();
        System.out.println();
        //        System.out.println(coniferRegionList);
    }

    public static void statisticSeqcnv(String path) {
        int sum = 0;
        int regionsInSimulated = 0;

        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String str = null;
            br.readLine();
            while ((str = br.readLine()) != null) {
                String[] features = str.split("\\s+");
                sum++;
                Region region = new Region(features[0], Integer.parseInt(features[1]),
                    Integer.parseInt(features[2]));
                if (judgeInSimulatedRegions(region)) {
                    seqcnvRegionList.add(region);
                    regionsInSimulated++;
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("SeqCNV count:" + sum + ", regionsInSimulated:" + regionsInSimulated);
        Collections.sort(seqcnvRegionList);

        for (Region region : simulatedRegionList) {
            seqcnvRegionMap.put(region, 0);
        }
        for (Region region : seqcnvRegionList) {
            Region ref = null;
            for (Map.Entry<Region, Integer> entry : seqcnvRegionMap.entrySet()) {
                if (entry.getKey().intersact(region)) {
                    ref = entry.getKey();
                    break;
                }
            }
            if (ref != null) {
                seqcnvRegionMap.put(ref, seqcnvRegionMap.get(ref) + 1);
            }
        }
        for (Map.Entry<Region, Integer> entry : seqcnvRegionMap.entrySet()) {
            System.out.println(entry.getKey() + "   " + entry.getValue());
        }
        System.out.println();
        System.out.println();
        System.out.println();
        //        System.out.println(seqcnvRegionList);
    }

    public static Set<Region> intersect(List<Region> region1, List<Region> region2) {
        Set<Region> result = new HashSet<Region>();
        for (Region regionItem1 : region1) {
            for (Region regionItem2 : region2) {
                int com = regionItem1.getChr().compareTo(regionItem2.getChr());
                if (com != 0)
                    continue;
                if (regionItem1.intersact(regionItem2))
                    result.add(regionItem1);
            }
        }
        return result;
    }
}
