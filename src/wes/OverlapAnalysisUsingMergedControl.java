package wes;

/**
 * Merge some bam files to regard as control, randomly select other samples to regard as cases.
 * 
 * <p>This overlap analysis use the strategy that if there is at least one bp overlaped between 
 * result and validation dataset(here we use CNV predicted results by conrad), we believe that 
 * the CNV event is predicted by the method.</p>
 * 
 * @author yorkchen
 * @since 2016-08-24
 *
 */
public class OverlapAnalysisUsingMergedControl {
	
	private String[] sampleCNVFileNames = {"NA10847_cnv.txt", "NA18967_cnv.txt", 
			"NA18973_cnv.txt", "NA18981_cnv.txt", "NA19131_cnv.txt"};
	
	private String conradCNVFilePaths = "E:\\WES";
	
	private String cnvResultFilePaths = "E:\\WES\\sample_merged_control_results";

	public static void main(String[] args) {

	}

}
