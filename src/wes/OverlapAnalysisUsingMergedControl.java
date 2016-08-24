package wes;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Merge some bam files to regard as control, randomly select other samples to regard as cases.
 * 
 * <p>
 * This overlap analysis use the strategy that if there is at least one bp overlaped between result
 * and validation dataset(here we use CNV predicted results by conrad), we believe that the CNV
 * event is predicted by the method.
 * </p>
 * 
 * @author yorkchen
 * @since 2016-08-24
 *
 */
public class OverlapAnalysisUsingMergedControl {

    private String[] sampleCNVFileNames = { "NA10847_cnv.txt", "NA18967_cnv.txt", "NA18973_cnv.txt",
                                            "NA18981_cnv.txt", "NA19131_cnv.txt" };

    private Path     conradCNVFilePaths = Paths.get("E:\\WES");

    private Path     cnvResultFilePaths = Paths.get("E:\\WES\\sample_merged_control_results");

    public static void main(String[] args) {

    }

    public void overlapAnalysis() throws IOException {
        List<Path> mergedControlPaths = Files.list(cnvResultFilePaths)
            .filter(path -> Files.isDirectory(path)).collect(Collectors.toList());
        for (Path mergedControlPath : mergedControlPaths) {
            String mergedControlName = mergedControlPath.getFileName().toString();
            List<Path> caseSampleResultPaths = Files.list(mergedControlPath)
                .filter(path -> Files.isDirectory(mergedControlPath)).collect(Collectors.toList());
            for (Path caseSampleResultPath : caseSampleResultPaths) {
                String caseSampleName = caseSampleResultPath.getFileName().toString();
                Path caseResultFilePath = caseSampleResultPath.resolve("report")
                    .resolve("CNV_report.txt");

            }

        }
    }

}
