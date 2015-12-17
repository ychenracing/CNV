package excavator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import utils.Region;

public class CopyTargetNotation {

    public static void main(String[] args) {
        copyTargetNotation("C:\\Users\\Administrator\\Desktop\\chr1_chr4.bed",
            "C:\\Users\\Administrator\\Desktop\\chr1_chr4_merged_3col.bed",
            "C:\\Users\\Administrator\\Desktop\\chr1_chr4_merged_notated.bed");
    }

    /**
     * copy target notation(columns start from 4 in source bed file) to target bed file.
     * @param srcBedPath
     * @param targetBedPath
     * @param outputBedPath
     */
    public static void copyTargetNotation(String srcBedPath, String targetBedPath,
                                          String outputBedPath) {
        File outputFile = new File(outputBedPath);
        if (!outputFile.exists()) {
            try {
                outputFile.createNewFile();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        List<String> srcLines = readLines(srcBedPath);
        List<String> targetLines = readLines(targetBedPath);
        try (BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile))) {
            int targetIndex = 0, srcIndex = 0;
            String srcLine = srcLines.get(srcIndex++);
            String[] srcFeature = srcLine.split("\\s+");
            String lastAppendString = "\t" + srcFeature[3] + "\t" + srcFeature[4] + "\t"
                                      + srcFeature[5] + "\n";
            String appendString = null;
            for (String targetLine : targetLines) {
                String[] targetFeature = targetLine.split("\\s+");
                Region targetRegion = new Region(targetFeature[0], targetFeature[1],
                    targetFeature[2]);
                while (srcIndex < srcLines.size()) {
                    srcLine = srcLines.get(srcIndex++);
                    srcFeature = srcLine.split("\\s+");
                    appendString = "\t" + srcFeature[3] + "\t" + srcFeature[4] + "\t"
                                   + srcFeature[5] + "\n";
                    Region srcRegion = new Region(srcFeature[0], srcFeature[1], srcFeature[2]);
                    if (!targetRegion.intersact(srcRegion)) {
                        String writeString = targetRegion + lastAppendString;
                        lastAppendString = appendString;
                        bw.write(writeString);
                        targetIndex++;
                        break;
                    } else {
                        srcIndex++;
                    }
                }
                if (srcIndex == srcLines.size()) {
                    String writeString = targetRegion + appendString;
                    bw.write(writeString);
                    targetIndex++;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Done!");
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
