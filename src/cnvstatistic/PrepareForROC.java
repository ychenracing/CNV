package cnvstatistic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.atomic.AtomicLong;

import utils.Pair;
import utils.Region;

/**
 * 之前计算FPR时，不知道TN，计算错误。尝试使用cnv region的长度(position/bp length)来计算各个工具的FPR。
 * @author racing
 * @version $Id: PrepareForROC.java, v 0.1 Nov 2, 2015 9:38:23 PM racing Exp $
 */
public class PrepareForROC {

    /** 所有reads总长度。应该使用bam文件里的reads总长度，但是怎么知道bam中reads的总长度？bed文件得不出reads总长度 */
    private long                allReadsLength     = 0;

    /** SeqCNV结果中case的chromosome文件夹，用以计算reads覆盖到的总长度 */
    private static final String CASE_CHR_FOLDER    = "/Users/racing/Lab/CNV/SeqCNV/ychen/cnvoutput3/tumor";
    /** SeqCNV结果中control的chromosome文件夹，用以计算reads覆盖到的总长度 */
    private static final String CONTROL_CHR_FOLDER = "/Users/racing/Lab/CNV/SeqCNV/ychen/cnvoutput3/control";

    /** 模拟的CNV区间的总长度 */
    private int                 simulatedLength    = 0;

    private Set<Region>         simulatedRegions   = new HashSet<>();
    private Set<Region>         seqcnvRegions      = new HashSet<>();
    private Set<Region>         coniferRegions     = new HashSet<>();
    private Set<Region>         cnvnatorRegions    = new HashSet<>();
    private Set<Region>         cnverRegions       = new HashSet<>();

    public static void main(String[] args) {
        PrepareForROC roc = new PrepareForROC();
        roc.calculateAllReadsLength();
        roc.readSimulatedLength("/Users/racing/Lab/(2015-09-22)CNV/simulatedRegions.txt");
        roc.calculateSeqCNVFPR("/Users/racing/Lab/(2015-09-22)CNV/SeqCNV/CNV_report.txt");
        roc.calculateConiferFPR("/Users/racing/Lab/(2015-09-22)CNV/Conifer/calls.txt");
        roc.calculateCNVnatorFPR(
            "/Users/racing/Lab/(2015-09-22)CNV/CNVnator/placenta_BAC_predict.txt");
        roc.calculateCNVerFPR("/Users/racing/Lab/(2015-09-22)CNV/CNVer");
    }

    /**
     * 从SeqCNV结果中tumor、control文件夹中读取所有reads的长度总和到变量allReadsLength中。
     */
    public void calculateAllReadsLength() {
        List<File> caseChrFileList = new ArrayList<>();
        caseChrFileList.addAll(Arrays.asList(new File(CASE_CHR_FOLDER).listFiles()));
        Set<String> chrFileNames = new HashSet<>();
        for (File item : caseChrFileList) {
            chrFileNames.add(item.getName());
        }

        AtomicLong sumLength = new AtomicLong(0);
        caseChrFileList.stream().forEach(chrFile -> {
            Pair<String, String> caseStartEnd = readFirstAndLastLine(chrFile);
            if (caseStartEnd == null)
                throw new RuntimeException("read first and last line error!");
            int caseStart = Math.abs(Integer.parseInt(caseStartEnd.getFirst()));
            int caseEnd = Math.abs(Integer.parseInt(caseStartEnd.getSecond()));

            Pair<String, String> controlStartEnd = readFirstAndLastLine(
                new File(CONTROL_CHR_FOLDER + File.separator + chrFile.getName()));
            if (controlStartEnd == null)
                throw new RuntimeException("read first and last line error!");
            int controlStart = Math.abs(Integer.parseInt(controlStartEnd.getFirst()));
            int controlEnd = Math.abs(Integer.parseInt(controlStartEnd.getSecond()));

            int chrStart = caseStart < controlStart ? caseStart : controlStart;
            int chrEnd = caseEnd > controlEnd ? caseEnd : controlEnd;
            sumLength.addAndGet(chrEnd - chrStart);
        });
        this.allReadsLength = sumLength.get();
    }

    /**
     * 读取simulate的区间的总长度到变量simulatedLength中。
     * @param simulatedRegionFilePath
     */
    public void readSimulatedLength(String simulatedRegionFilePath) {
        int sumLength = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(simulatedRegionFilePath))) {
            String input = null;
            while ((input = br.readLine()) != null) {
                String[] feature = input.split("\\s+");
                sumLength += Integer.parseInt(feature[2]) - Integer.parseInt(feature[1]);
                this.simulatedRegions.add(new Region(feature[0], Integer.parseInt(feature[1]),
                    Integer.parseInt(feature[2])));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        this.simulatedLength = sumLength;
    }

    /**
     * 根据输入的SeqCNV预测结果文件的路径，根据region的长度来计算FPR等指标。
     * @param seqcnvResultFilePath
     */
    public void calculateSeqCNVFPR(String seqcnvResultFilePath) {
        long tnLength = 0;
        long fpLength = 0;
        long fnLength = 0;
        long tpLength = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(seqcnvResultFilePath))) {
            String input = br.readLine(); // discard header
            while ((input = br.readLine()) != null) {
                String[] feature = input.split("\\s+");
                Region predictedRegion = new Region(feature[0], feature[1], feature[2]);
                this.seqcnvRegions.add(predictedRegion);
                boolean isIntersacted = false;
                for (Region simulatedRegion : simulatedRegions) {
                    if (predictedRegion.intersact(simulatedRegion)) {
                        isIntersacted = true;
                        // predictedRegion完全包含在simulatedRegion中
                        if (predictedRegion.getStart() >= simulatedRegion.getStart()
                            && predictedRegion.getEnd() <= simulatedRegion.getEnd()) {
                            tpLength += predictedRegion.getEnd() - predictedRegion.getStart();
                        }
                        // predictedRegion的右边与simulatedRegion左边重叠
                        else if (predictedRegion.getStart() < simulatedRegion.getStart()
                                 && predictedRegion.getEnd() <= simulatedRegion.getEnd()) {
                            tpLength += predictedRegion.getEnd() - simulatedRegion.getStart();
                            fpLength += simulatedRegion.getStart() - predictedRegion.getStart();
                        }
                        // predictedRegion的左边与simulatedRegion右边重叠
                        else if (predictedRegion.getStart() >= simulatedRegion.getStart()
                                 && predictedRegion.getEnd() >= simulatedRegion.getEnd()) {
                            tpLength += simulatedRegion.getEnd() - predictedRegion.getStart();
                            fpLength += predictedRegion.getEnd() - simulatedRegion.getEnd();
                        }
                        // simulatedRegion完全包含在predictedRegion中
                        else if (predictedRegion.getStart() <= simulatedRegion.getStart()
                                 && predictedRegion.getEnd() >= simulatedRegion.getEnd()) {
                            tpLength += simulatedRegion.getEnd() - simulatedRegion.getStart();
                            fpLength += simulatedRegion.getStart() - predictedRegion.getStart();
                            fpLength += predictedRegion.getEnd() - simulatedRegion.getEnd();
                        }
                        break;
                    }
                }

                if (!isIntersacted) {
                    fpLength += predictedRegion.getEnd() - predictedRegion.getStart();
                }
            }

            fnLength = this.simulatedLength - tpLength;
            tnLength = this.allReadsLength - this.simulatedLength - fpLength;

            System.out.println("seqcnvPredictedLength: " + this.seqcnvRegions.stream()
                .map(region -> region.getEnd() - region.getStart()).reduce(0L, (x, y) -> x + y));
            System.out.println("allReadsLength: " + this.allReadsLength);
            System.out.println("simulatedLength: " + this.simulatedLength);
            System.out.println("fpLength: " + fpLength);
            System.out.println("tnLength: " + tnLength);
            System.out.println("fnLength: " + fnLength);
            System.out.println("tpLength: " + tpLength);

            // FPR = FP / (FP + TN)
            float fpr = (float) fpLength / (fpLength + tnLength);
            // Sensitivity = Recall = TP / (TP + FN)
            float recall = (float) tpLength / (tpLength + fnLength);
            // Specificity = 1 - FPR
            float specificity = 1 - fpr;

            System.out.println("SeqCNV FPR = FP / (FP + TN) = " + fpr);
            System.out.println("SeqCNV Sensitivity = TP / (TP + FN) = " + recall);
            System.out.println("SeqCNV Specificity = 1 - FPR = " + specificity);
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private Pair<String, String> readFirstAndLastLine(File file) {
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String input = br.readLine();
            String start = input;
            String lastLine = null;
            while ((input = br.readLine()) != null) {
                lastLine = input;
            }
            String end = lastLine;
            return new Pair<String, String>(start, end);
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    /**
     * 输入CoNIFER预测结果文件路径，根据length计算FPR
     * @param coniferResultFilePath
     */
    public void calculateConiferFPR(String coniferResultFilePath) {
        long tnLength = 0;
        long fpLength = 0;
        long fnLength = 0;
        long tpLength = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(coniferResultFilePath))) {
            String line = br.readLine(); // discard header
            while ((line = br.readLine()) != null) {
                if (line.startsWith("placenta_BAC")) {
                    String[] feature = line.split("\\s+");
                    Region predictedRegion = new Region(feature[1], feature[2], feature[3]);
                    this.coniferRegions.add(predictedRegion);
                    boolean isIntersacted = false;
                    for (Region simulatedRegion : simulatedRegions) {
                        if (predictedRegion.intersact(simulatedRegion)) {
                            isIntersacted = true;
                            // predictedRegion完全包含在simulatedRegion中
                            if (predictedRegion.getStart() >= simulatedRegion.getStart()
                                && predictedRegion.getEnd() <= simulatedRegion.getEnd()) {
                                tpLength += predictedRegion.getEnd() - predictedRegion.getStart();
                            }
                            // predictedRegion的右边与simulatedRegion左边重叠
                            else if (predictedRegion.getStart() < simulatedRegion.getStart()
                                     && predictedRegion.getEnd() <= simulatedRegion.getEnd()) {
                                tpLength += predictedRegion.getEnd() - simulatedRegion.getStart();
                                fpLength += simulatedRegion.getStart() - predictedRegion.getStart();
                            }
                            // predictedRegion的左边与simulatedRegion右边重叠
                            else if (predictedRegion.getStart() >= simulatedRegion.getStart()
                                     && predictedRegion.getEnd() >= simulatedRegion.getEnd()) {
                                tpLength += simulatedRegion.getEnd() - predictedRegion.getStart();
                                fpLength += predictedRegion.getEnd() - simulatedRegion.getEnd();
                            }
                            // simulatedRegion完全包含在predictedRegion中
                            else if (predictedRegion.getStart() <= simulatedRegion.getStart()
                                     && predictedRegion.getEnd() >= simulatedRegion.getEnd()) {
                                tpLength += simulatedRegion.getEnd() - simulatedRegion.getStart();
                                fpLength += simulatedRegion.getStart() - predictedRegion.getStart();
                                fpLength += predictedRegion.getEnd() - simulatedRegion.getEnd();
                            }
                            break;
                        }
                    }

                    if (!isIntersacted) {
                        fpLength += predictedRegion.getEnd() - predictedRegion.getStart();
                    }
                }
            }

            fnLength = this.simulatedLength - tpLength;
            tnLength = this.allReadsLength - this.simulatedLength - fpLength;

            System.out.println("coniferPredictedLength: " + this.coniferRegions.stream()
                .map(region -> region.getEnd() - region.getStart()).reduce(0L, (x, y) -> x + y));
            System.out.println("allReadsLength: " + this.allReadsLength);
            System.out.println("simulatedLength: " + this.simulatedLength);
            System.out.println("fpLength: " + fpLength);
            System.out.println("tnLength: " + tnLength);
            System.out.println("fnLength: " + fnLength);
            System.out.println("tpLength: " + tpLength);

            // FPR = FP / (FP + TN)
            float fpr = (float) fpLength / (fpLength + tnLength);
            // Sensitivity = Recall = TP / (TP + FN)
            float recall = (float) tpLength / (tpLength + fnLength);
            // Specificity = 1 - FPR
            float specificity = 1 - fpr;

            System.out.println("CoNIFER FPR = FP / (FP + TN) = " + fpr);
            System.out.println("CoNIFER Sensitivity = TP / (TP + FN) = " + recall);
            System.out.println("CoNIFER Specificity = 1 - FPR = " + specificity);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 输入CNVnator预测结果文件路径，根据length计算FPR
     * @param cnvnatorResultFilePath
     */
    public void calculateCNVnatorFPR(String cnvnatorResultFilePath) {
        long tnLength = 0;
        long fpLength = 0;
        long fnLength = 0;
        long tpLength = 0;
        try (BufferedReader br = new BufferedReader(new FileReader(cnvnatorResultFilePath))) {
            String line = null; // discard header
            while ((line = br.readLine()) != null) {
                String[] feature = line.split("\\s+")[1].replace("-", ":").split(":");
                Region predictedRegion = new Region("chr" + feature[0], feature[1], feature[2]);
                this.cnvnatorRegions.add(predictedRegion);
                boolean isIntersacted = false;
                for (Region simulatedRegion : simulatedRegions) {
                    if (predictedRegion.intersact(simulatedRegion)) {
                        isIntersacted = true;
                        // predictedRegion完全包含在simulatedRegion中
                        if (predictedRegion.getStart() >= simulatedRegion.getStart()
                            && predictedRegion.getEnd() <= simulatedRegion.getEnd()) {
                            tpLength += predictedRegion.getEnd() - predictedRegion.getStart();
                        }
                        // predictedRegion的右边与simulatedRegion左边重叠
                        else if (predictedRegion.getStart() < simulatedRegion.getStart()
                                 && predictedRegion.getEnd() <= simulatedRegion.getEnd()) {
                            tpLength += predictedRegion.getEnd() - simulatedRegion.getStart();
                            fpLength += simulatedRegion.getStart() - predictedRegion.getStart();
                        }
                        // predictedRegion的左边与simulatedRegion右边重叠
                        else if (predictedRegion.getStart() >= simulatedRegion.getStart()
                                 && predictedRegion.getEnd() >= simulatedRegion.getEnd()) {
                            tpLength += simulatedRegion.getEnd() - predictedRegion.getStart();
                            fpLength += predictedRegion.getEnd() - simulatedRegion.getEnd();
                        }
                        // simulatedRegion完全包含在predictedRegion中
                        else if (predictedRegion.getStart() <= simulatedRegion.getStart()
                                 && predictedRegion.getEnd() >= simulatedRegion.getEnd()) {
                            tpLength += simulatedRegion.getEnd() - simulatedRegion.getStart();
                            fpLength += simulatedRegion.getStart() - predictedRegion.getStart();
                            fpLength += predictedRegion.getEnd() - simulatedRegion.getEnd();
                        }
                        break;
                    }
                }

                if (!isIntersacted) {
                    fpLength += predictedRegion.getEnd() - predictedRegion.getStart();
                }
            }

            fnLength = this.simulatedLength - tpLength;
            tnLength = this.allReadsLength - this.simulatedLength - fpLength;

            System.out.println("cnvnatorPredictedLength: " + this.cnvnatorRegions.stream()
                .map(region -> region.getEnd() - region.getStart()).reduce(0L, (x, y) -> x + y));
            System.out.println("allReadsLength: " + this.allReadsLength);
            System.out.println("simulatedLength: " + this.simulatedLength);
            System.out.println("fpLength: " + fpLength);
            System.out.println("tnLength: " + tnLength);
            System.out.println("fnLength: " + fnLength);
            System.out.println("tpLength: " + tpLength);

            // FPR = FP / (FP + TN)
            float fpr = (float) fpLength / (fpLength + tnLength);
            // Sensitivity = Recall = TP / (TP + FN)
            float recall = (float) tpLength / (tpLength + fnLength);
            // Specificity = 1 - FPR
            float specificity = 1 - fpr;

            System.out.println("CNVnator FPR = FP / (FP + TN) = " + fpr);
            System.out.println("CNVnator Sensitivity = TP / (TP + FN) = " + recall);
            System.out.println("CNVnator Specificity = 1 - FPR = " + specificity);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 输入CNVnator预测结果目录路径，根据length计算FPR
     * @param cnverResultFolderPath
     */
    public void calculateCNVerFPR(String cnverResultFolderPath) {
        AtomicLong tnLength = new AtomicLong(0);
        AtomicLong fpLength = new AtomicLong(0);
        AtomicLong fnLength = new AtomicLong(0);
        AtomicLong tpLength = new AtomicLong(0);
        File[] resultFiles = new File(cnverResultFolderPath).listFiles();
        List<File> list = new ArrayList<>();
        list.addAll(Arrays.asList(resultFiles));
        list.stream().forEach(resultFileItem -> {
            try (BufferedReader br = new BufferedReader(new FileReader(resultFileItem))) {
                String line = null; // discard header
                while ((line = br.readLine()) != null) {
                    String[] feature = line.split("\\s+");
                    Region predictedRegion = new Region(feature[0], feature[1], feature[2]);
                    this.cnverRegions.add(predictedRegion);
                    boolean isIntersacted = false;
                    for (Region simulatedRegion : simulatedRegions) {
                        if (predictedRegion.intersact(simulatedRegion)) {
                            isIntersacted = true;
                            // predictedRegion完全包含在simulatedRegion中
                            if (predictedRegion.getStart() >= simulatedRegion.getStart()
                                && predictedRegion.getEnd() <= simulatedRegion.getEnd()) {
                                tpLength.getAndAdd(
                                    predictedRegion.getEnd() - predictedRegion.getStart());
                            }
                            // predictedRegion的右边与simulatedRegion左边重叠
                            else if (predictedRegion.getStart() < simulatedRegion.getStart()
                                     && predictedRegion.getEnd() <= simulatedRegion.getEnd()) {
                                tpLength.getAndAdd(
                                    predictedRegion.getEnd() - simulatedRegion.getStart());
                                fpLength.getAndAdd(
                                    simulatedRegion.getStart() - predictedRegion.getStart());
                            }
                            // predictedRegion的左边与simulatedRegion右边重叠
                            else if (predictedRegion.getStart() >= simulatedRegion.getStart()
                                     && predictedRegion.getEnd() >= simulatedRegion.getEnd()) {
                                tpLength.getAndAdd(
                                    simulatedRegion.getEnd() - predictedRegion.getStart());
                                fpLength
                                    .getAndAdd(predictedRegion.getEnd() - simulatedRegion.getEnd());
                            }
                            // simulatedRegion完全包含在predictedRegion中
                            else if (predictedRegion.getStart() <= simulatedRegion.getStart()
                                     && predictedRegion.getEnd() >= simulatedRegion.getEnd()) {
                                tpLength.getAndAdd(
                                    simulatedRegion.getEnd() - simulatedRegion.getStart());
                                fpLength.getAndAdd(
                                    simulatedRegion.getStart() - predictedRegion.getStart());
                                fpLength
                                    .getAndAdd(predictedRegion.getEnd() - simulatedRegion.getEnd());
                            }
                            break;
                        }
                    }

                    if (!isIntersacted) {
                        fpLength.getAndAdd(predictedRegion.getEnd() - predictedRegion.getStart());
                    }
                }

            } catch (IOException e) {
                e.printStackTrace();
            }
        });
        fnLength.set(this.simulatedLength - tpLength.get());
        tnLength.set(this.allReadsLength - this.simulatedLength - fpLength.get());

        System.out.println("cnverPredictedLength: " + this.cnverRegions.stream()
            .map(region -> region.getEnd() - region.getStart()).reduce(0L, (x, y) -> x + y));
        System.out.println("allReadsLength: " + this.allReadsLength);
        System.out.println("simulatedLength: " + this.simulatedLength);
        System.out.println("fpLength: " + fpLength);
        System.out.println("tnLength: " + tnLength);
        System.out.println("fnLength: " + fnLength);
        System.out.println("tpLength: " + tpLength);

        // FPR = FP / (FP + TN)
        float fpr = (float) fpLength.get() / (fpLength.get() + tnLength.get());
        // Sensitivity = Recall = TP / (TP + FN)
        float recall = (float) tpLength.get() / (tpLength.get() + fnLength.get());
        // Specificity = 1 - FPR
        float specificity = 1 - fpr;

        System.out.println("CNVer FPR = FP / (FP + TN) = " + fpr);
        System.out.println("CNVer Sensitivity = TP / (TP + FN) = " + recall);
        System.out.println("CNVer Specificity = 1 - FPR = " + specificity);
    }
}
