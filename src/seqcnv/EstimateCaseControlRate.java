package seqcnv;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.stream.Stream;

import utils.Pair;
import utils.Region;

public class EstimateCaseControlRate {

    private static List<Region>                                     simulatedRegionList = new ArrayList<Region>();

    private static Map<String, String>                              tumorChrStartEnd    = new HashMap<String, String>();
    private static Map<String, String>                              controlChrStartEnd  = new HashMap<String, String>();

    /** 每条染色体中case的reads数量和control的reads数量 */
    private static Map<String, Pair<Integer, Integer>>              chrReadNum          = new HashMap<String, Pair<Integer, Integer>>();
    /** 每个region中case的reads数量和control的reads数量 */
    private static Map<Region, Pair<Float, Pair<Integer, Integer>>> regionReadNum       = new LinkedHashMap<Region, Pair<Float, Pair<Integer, Integer>>>();

    public static void main(String[] args) {
        recordStartEndPoint("/Users/racing/Lab/CNV/SeqCNV/ychen/cnvoutput3/tumor", true);
        recordStartEndPoint("/Users/racing/Lab/CNV/SeqCNV/ychen/cnvoutput3/control", false);
        countReadsForChrs();
        readSimulatedRegions("/Users/racing/Lab/(2015-09-22)CNV/simulatedRegions.txt");
        countRegionReadsRatio("/Users/racing/Lab/CNV/SeqCNV_report.txt");
        recordCaseControlReadsRatio("/Users/racing/Lab/CNV");
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

    private static void readSimulatedRegions(String path) {
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String str = null;
            while ((str = br.readLine()) != null) {
                String[] features = str.split("\\s+");
                simulatedRegionList.add(new Region(features[0], Integer.parseInt(features[1]),
                    Integer.parseInt(features[2])));
            }
            System.out.println("simulated regions are read into simulatedRegionList!");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 把每条染色体的reads数量当做library size计算case-control reads ratio，输出到chrLibrarySize.txt中。
     * 把所有染色体的reads数量当做library size计算case-control reads ratio，输出到totalLibrarySize.txt中。
     * @param folder
     */
    public static void recordCaseControlReadsRatio(String folder) {
        File chrLibrarySize = new File(folder + File.separator + "chrLibrarySize.txt");
        File totalLibrarySize = new File(folder + File.separator + "totalLibrarySize.txt");
        try {
            if (!chrLibrarySize.exists())
                chrLibrarySize.createNewFile();
            if (!totalLibrarySize.exists())
                totalLibrarySize.createNewFile();
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        // 首先计算每条染色体的reads数量为library size的正则化
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(chrLibrarySize))) {
            String header = "#chr\tstart\tend\tcase_reads_num\tcontrol_reads_num\tTPorFP\tpre_ratio\tratio\n";
            bw.write(header, 0, header.length());
            System.out.println(header);
            regionReadNum.entrySet().stream().forEach(entry -> {
                Region region = entry.getKey();
                String trueOrFalsePositive = "FP";
                for (Region item : simulatedRegionList) {
                    if (region.intersact(item)) {
                        trueOrFalsePositive = "TP";
                        break;
                    }
                }
                StringBuilder writeString = new StringBuilder(region.toString()).append("\t");
                Pair<Float, Pair<Integer, Integer>> pair = entry.getValue();
                Pair<Integer, Integer> readsNum = pair.getSecond();

                int caseReadsNum = readsNum.getFirst();
                int controlReadsNum = readsNum.getSecond();
                int chrCaseReadsNum = chrReadNum.get(entry.getKey().getChr()).getFirst();
                int chrControlReadsNum = chrReadNum.get(entry.getKey().getChr()).getSecond();
                float seqcnvPredictRatio = pair.getFirst();

                float readsRatio = (float) caseReadsNum / controlReadsNum;
                float librarySizeRatio = (float) chrControlReadsNum / chrCaseReadsNum;
                float ratio = readsRatio * librarySizeRatio;

                System.out.println(region + "\tcaseReadsNum:" + caseReadsNum + "\tcontrolReadsNum:"
                                   + controlReadsNum);
                System.out.println("chrCaseReadsSum:" + chrCaseReadsNum + "\tchrControlReadsSum:"
                                   + chrControlReadsNum);
                System.out.println("chrLibrarySizeRatio:" + librarySizeRatio);

                writeString.append(caseReadsNum).append("\t");
                writeString.append(controlReadsNum).append("\t");
                writeString.append(trueOrFalsePositive).append("\t");
                writeString.append(seqcnvPredictRatio).append("\t");
                writeString.append(String.format("%.3f", ratio)).append("\n");
                try {
                    bw.write(writeString.toString(), 0, writeString.length());
                    System.out.println(writeString.toString().trim());
                } catch (Exception e) {
                    e.printStackTrace();
                }
            });
        } catch (IOException ex) {
            ex.printStackTrace();
        }

        System.out.println();
        System.out.println();

        // 再计算所有染色体的reads数量为library size的正则化
        int caseReadsSum = chrReadNum.entrySet().stream().map(entry -> entry.getValue().getFirst())
            .reduce(0, (x, y) -> x + y);
        int controlReadsSum = chrReadNum.entrySet().stream()
            .map(entry -> entry.getValue().getSecond()).reduce(0, (x, y) -> x + y);
        float librarySizeRatio = (float) controlReadsSum / caseReadsSum;

        System.out.println("caseReadsSum:" + caseReadsSum + "\tcontrolReadsSum:" + controlReadsSum);
        System.out.println("librarySizeRatio:" + librarySizeRatio);

        try (BufferedWriter bw = new BufferedWriter(new FileWriter(totalLibrarySize))) {
            String header = "#chr\tstart\tend\tcase_reads_num\tcontrol_reads_num\tTPorFP\tpre_ratio\tratio\n";
            bw.write(header, 0, header.length());
            System.out.println(header);
            regionReadNum.entrySet().stream().forEach(entry -> {
                Region region = entry.getKey();
                String trueOrFalsePositive = "FP";
                for (Region item : simulatedRegionList) {
                    if (region.intersact(item)) {
                        trueOrFalsePositive = "TP";
                        break;
                    }
                }
                StringBuilder writeString = new StringBuilder(region.toString()).append("\t");

                Pair<Float, Pair<Integer, Integer>> pair = entry.getValue();
                Pair<Integer, Integer> readsNum = pair.getSecond();
                int caseReadsNum = readsNum.getFirst();
                int controlReadsNum = readsNum.getSecond();
                float readsRatio = (float) caseReadsNum / controlReadsNum;
                float ratio = readsRatio * librarySizeRatio;
                float seqcnvPredictRatio = pair.getFirst();

                writeString.append(caseReadsNum).append("\t");
                writeString.append(controlReadsNum).append("\t");
                writeString.append(seqcnvPredictRatio).append("\t");
                writeString.append(trueOrFalsePositive).append("\t");
                writeString.append(String.format("%.3f", ratio)).append("\n");
                try {
                    bw.write(writeString.toString(), 0, writeString.length());
                    System.out.println(writeString.toString().trim());
                } catch (Exception e) {
                    e.printStackTrace();
                }
            });
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        System.out.println();
        System.out.println();
    }

    /**
     * 计算SeqCNV的输出中，每个假阳的区间的case的reads数量和control的reads数量，保存在regionReadNum中。
     * @param cnvResultFilePath
     */
    public static void countRegionReadsRatio(String cnvResultFilePath) {
        File file = new File(cnvResultFilePath);
        try (Scanner scanner = new Scanner(new FileInputStream(file))) {
            String input = null;
            while (scanner.hasNextLine()) {
                input = scanner.nextLine();
                String[] feature = input.split("\\s+");
                Region region = new Region(feature[0], Integer.parseInt(feature[1]),
                    Integer.parseInt(feature[2]));

                String regionString = feature[0] + ":" + numberConvertToString(feature[1]) + "-"
                                      + numberConvertToString(feature[2]);

                int controlReadNum = 0;
                int caseReadNum = 0;
                try {
                    String[] controlCmds = { "/bin/bash", "-cl",
                                             "samtools view /Users/racing/Lab/CNV/CNVnator/ychen/placenta.local_realigned.sorted.bam "
                                                                 + regionString + " | wc -l" };
                    Process controlProcess = new ProcessBuilder(controlCmds).start();
                    try (BufferedReader br = new BufferedReader(
                        new InputStreamReader(controlProcess.getInputStream()))) {
                        controlReadNum = Integer.parseInt(br.readLine().trim());
                    }
                    String[] caseCmds = { "/bin/bash", "-cl",
                                          "samtools view /Users/racing/Lab/CNV/CNVnator/ychen/placenta_BAC.local_realigned.sorted.bam "
                                                              + regionString + " | wc -l" };
                    Process caseProcess = Runtime.getRuntime().exec(caseCmds);
                    try (BufferedReader br = new BufferedReader(
                        new InputStreamReader(caseProcess.getInputStream()))) {
                        caseReadNum = Integer.parseInt(br.readLine().trim());
                    }

                    System.out.println("counting region reads ratio for " + region + " finished!");
                    regionReadNum.put(region,
                        new Pair<Float, Pair<Integer, Integer>>(Float.parseFloat(feature[6]),
                            new Pair<Integer, Integer>(controlReadNum, caseReadNum)));
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 计算每条染色体上的case中的reads数量和control中的reads数量
     */
    public static void countReadsForChrs() {
        tumorChrStartEnd.entrySet().stream().forEach(entry -> {
            int controlReadNum = 0;
            int caseReadNum = 0;
            try {
                String[] controlCmds = { "/bin/bash", "-cl",
                                         "samtools view /Users/racing/Lab/CNV/CNVnator/ychen/placenta.local_realigned.sorted.bam "
                                                             + controlChrStartEnd.get(
                                                                 entry.getKey())
                                                             + " | wc -l" };
                Process controlProcess = new ProcessBuilder(controlCmds).start();
                try (BufferedReader br = new BufferedReader(
                    new InputStreamReader(controlProcess.getInputStream()))) {
                    String input = null;
                    while ((input = br.readLine()) != null) {
                        controlReadNum = Integer.parseInt(input.trim());
                    }
                }
                String[] caseCmds = { "/bin/bash", "-cl",
                                      "samtools view /Users/racing/Lab/CNV/CNVnator/ychen/placenta_BAC.local_realigned.sorted.bam "
                                                          + entry.getValue() + " | wc -l" };
                Process caseProcess = Runtime.getRuntime().exec(caseCmds);
                try (BufferedReader br = new BufferedReader(
                    new InputStreamReader(caseProcess.getInputStream()))) {
                    String input = null;
                    while ((input = br.readLine()) != null) {
                        caseReadNum = Integer.parseInt(input.trim());
                    }
                }
                System.out.println(
                    "counting case control reads number for " + entry.getKey() + " finished!");
                chrReadNum.put(entry.getKey(),
                    new Pair<Integer, Integer>(caseReadNum, controlReadNum));
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
        System.out.println("countReadsForChrs finished!");
    }

    /**
     * 计算每条染色体的其实位置和终止位置，添加到tumorChrStartEnd和controlChrStartEnd中。
     * 其中，每一项的格式为：chr1:13,160-249,093,797
     * @param chrFolder
     * @param tumor
     */
    public static void recordStartEndPoint(String chrFolder, boolean tumor) {
        File folderFile = new File(chrFolder);
        File[] chrs = folderFile.listFiles();
        List<File> fileList = new ArrayList<File>(Arrays.asList(chrs));

        Stream<File> streams = fileList.stream();
        streams.filter(file -> {
            return file.getName().startsWith("chr");
        }).forEach(file -> {
            try (BufferedReader br = new BufferedReader(new FileReader(file))) {
                String input = numberConvertToString(br.readLine());
                String start = input;
                String lastLine = null;
                while ((input = br.readLine()) != null) {
                    lastLine = input;
                }
                String end = numberConvertToString(lastLine);
                if (tumor) {
                    tumorChrStartEnd.put(file.getName(), file.getName() + ":" + start + "-" + end);
                } else {
                    controlChrStartEnd.put(file.getName(),
                        file.getName() + ":" + start + "-" + end);
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        });
        System.out
            .println("recording " + (tumor ? "tumor" : "control") + " start end point finished!");
    }

    /**
     * 把字符串表示的数字转化为xxx,xxx,xxx格式的字符串，数字前可能有负号。
     * @param str
     * @return
     */
    private static String numberConvertToString(String str) {
        str = str.charAt(0) == '-' ? str.substring(1) : str;
        int commaNum = (int) Math.floor((double) (str.length() - 1) / 3);
        int remainder = str.length() % 3;
        StringBuilder result = new StringBuilder();
        int index = 0;
        if (remainder == 0) {
            for (int i = 0; i < 3; i++) {
                result.append(str.charAt(index++));
            }
        } else {
            while (index < remainder) {
                result.append(str.charAt(index++));
            }
        }
        for (int i = 0; i < commaNum; i++) {
            result.append(",");
            for (int j = 0; j < 3; j++) {
                result.append(str.charAt(index++));
            }
        }
        return result.toString();
    }
}
