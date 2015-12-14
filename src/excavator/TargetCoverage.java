package excavator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;

public class TargetCoverage {

    private static String placetaBamFilePath    = null;
    private static String placetaBacBamFilePath = null;

    public static void main(String[] args) {

        if (args.length < 4) {
            usage();
            return;
        }
        
        placetaBamFilePath = args[0];
        placetaBacBamFilePath = args[1];
        String bedFilePath = args[2];
        String outputFilePath = args[3];
        getCoverageFromTarget(bedFilePath, outputFilePath);
        System.out.println("done!");

    }

    /**
     * 计算bed中所有target region上的case和control的reads数量，输出到结果文件中。
     * @param bedFilePath
     * @param outputFilePath
     */
    public static void getCoverageFromTarget(String bedFilePath, String outputFilePath) {
        if(!new File(outputFilePath).exists()){
            try {
                new File(outputFilePath).createNewFile();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        try (BufferedReader br = new BufferedReader(new FileReader(bedFilePath));
                BufferedWriter bw = new BufferedWriter(new FileWriter(outputFilePath))) {
            bw.write("target\tplacenta_BAC\tplacenta\n");
            String readLine = null;
            while ((readLine = br.readLine()) != null) {
                String[] feature = readLine.split("\\s+");
                feature[0] = feature[0].length() < 3 ? "chr" + feature[0] : feature[0];
                String samtoolsRegion = feature[0] + ":" + numberConvertToString(feature[1]) + "-"
                                        + numberConvertToString(feature[2]);
                String samtoolsViewCmd = "samtools view " + placetaBamFilePath + " "
                                         + samtoolsRegion + " | wc -l";
                String[] execCmds = { "/bin/bash", "-cl", samtoolsViewCmd };
                Process subprocess = Runtime.getRuntime().exec(execCmds);
                BufferedReader subprocessBr = new BufferedReader(
                    new InputStreamReader(subprocess.getInputStream()));
                String placentaReadsCount = subprocessBr.readLine().trim();

                if (subprocess.waitFor() != 0) {
                    System.out
                        .println("exception region:" + samtoolsRegion + " for placenta" + "!!!!!!");
                } else {
                    System.out.println(samtoolsRegion + " for placenta finished");
                }
                subprocessBr.close();
                samtoolsViewCmd = "samtools view " + placetaBacBamFilePath + " " + samtoolsRegion
                                  + " | wc -l";
                execCmds[2] = samtoolsViewCmd;
                subprocess = Runtime.getRuntime().exec(execCmds);
                subprocessBr = new BufferedReader(
                    new InputStreamReader(subprocess.getInputStream()));
                String placentaBacReadsCount = subprocessBr.readLine().trim();
                if (subprocess.waitFor() != 0) {
                    System.out.println(
                        "exception region:" + samtoolsRegion + " for placenta_BAC" + "!!!!!!");
                } else {
                    System.out.println(samtoolsRegion + " for placenta_BAC finished");
                }
                subprocessBr.close();
                bw.write(feature[0] + ":" + feature[1] + "-" + feature[2] + "\t"
                         + placentaReadsCount + "\t" + placentaBacReadsCount);
            }
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
    }

    /**
     * 把字符串表示的数字转化为xxx,xxx,xxx格式的字符串，数字前可能有负号，如果有负号，负号会被剪切掉。
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

    public static void usage() {
        System.out.println(
            "Usage:\n\tjava TargetCoverage [caseBamFilePath] [controlBamFilePath] [bedFilePath] [resultFilePath]");
    }

}
