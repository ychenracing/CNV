package wes;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import utils.Region;
import utils.Region.CNVType;

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

    private static String[] sampleCNVFileNames = { "NA10847_cnv.txt", "NA18967_cnv.txt", 
    	"NA18973_cnv.txt", "NA18981_cnv.txt", "NA19131_cnv.txt" };

    private static Path     conradCNVFilePaths = Paths.get("E:\\WES");

    private static Path     cnvResultFilePaths = Paths.get("E:\\WES\\sample_merged_control_results");
    
    private static Map<String, List<Region>> conradCNVEventsMap = new HashMap<String, List<Region>>();
    
    private static Map<String, Map<String, List<String>>> formattedOutputResult = new LinkedHashMap<String, Map<String,List<String>>>();
    
    private static double seqcnvGainCutoff = 1.5;

    private static double seqcnvLossCutoff = 0.5;
    

    public static void main(String[] args) throws IOException {
//    	analyzeSeqCNVUsingConradCNVs();
    	
//    	Path resultPath = Paths.get("E:\\WES\\sample_merged_control_results\\NA12751_NA18999_NA18959_NA19153_NA12761\\NA10847\\report\\CNV_report.txt");
//    	Path filteredResultPath = Paths.get("E:\\WES\\sample_merged_control_results\\NA12751_NA18999_NA18959_NA19153_NA12761\\NA10847\\report\\filtered_CNV_report.txt");
//    	Path bedFilePath = Paths.get("C:\\Users\\yorkchen\\Desktop\\hglft_genome_3155_fe5330.bed");
//    	filterSeqCNVResultWithDetectingBed(resultPath, filteredResultPath, bedFilePath);
    	
    	Path filteredResultPath = Paths.get("E:\\WES\\sample_merged_control_results\\NA12751_NA18999_NA18959_NA19153_NA12761\\NA10847");
    	Path mergedControlFilePath = Paths.get("E:\\WES\\NA10847__NA12751_NA18999_NA18959_NA19153_NA12761.txt");
    	analyzeSeqCNVOneSampleUsingDetectingBedWithMergedControl(filteredResultPath, mergedControlFilePath);
    }

    public static void analyzeSeqCNVUsingConradCNVs() throws IOException {
        List<Path> mergedControlPaths = Files.list(cnvResultFilePaths)
            .filter(path -> Files.isDirectory(path)).collect(Collectors.toList());
        for (Path mergedControlPath : mergedControlPaths) {
            String mergedControlName = mergedControlPath.getFileName().toString();
            List<Path> caseSampleResultPaths = Files.list(mergedControlPath)
                .filter(path -> Files.isDirectory(path)).collect(Collectors.toList());
            for (Path caseSampleResultPath : caseSampleResultPaths) {
                String caseSampleName = caseSampleResultPath.getFileName().toString();
                Path caseResultFilePath = caseSampleResultPath.resolve("report")
                    .resolve("CNV_report.txt");
                List<String> sampleResultCNVEventsLines = Files.readAllLines(caseResultFilePath);
                sampleResultCNVEventsLines = sampleResultCNVEventsLines.stream().filter(line -> {
                	String[] cnvFeatures = line.split("\\s+");
                	long cnvStartPos = Long.parseLong(cnvFeatures[1]); 
                	long cnvEndPos = Long.parseLong(cnvFeatures[2]); 
                	Float cnvRatio = Float.parseFloat(cnvFeatures[cnvFeatures.length - 1]);
                	return (cnvRatio > seqcnvGainCutoff || cnvRatio < seqcnvLossCutoff) 
                			&& (cnvEndPos - cnvStartPos) > 1000;
                }).collect(Collectors.toList());
                if (!conradCNVEventsMap.containsKey(caseSampleName)) {
                    Optional<String> sampleCNVFileName = Stream.of(sampleCNVFileNames)
                    		.filter(fileName -> fileName.startsWith(caseSampleName)).findFirst();
                    Path sampleCNVFilePath = null;
                    if (sampleCNVFileName.isPresent()) {
                    	sampleCNVFilePath = conradCNVFilePaths.resolve(sampleCNVFileName.get());
                    }
                    List<String> sampleCNVEventsLines = Files.readAllLines(sampleCNVFilePath);
                    List<Region> sampleCNVEvents = new ArrayList<Region>();
                    sampleCNVEventsLines.remove(0);  // drop header
                    for (String sampleCNVEventLine : sampleCNVEventsLines) {
                    	String[] cnvFeatures = sampleCNVEventLine.split("\\s+");
                    	String chr = cnvFeatures[0].equals("1") ? "chr" + cnvFeatures[0] : 
                    		cnvFeatures[0];
                    	CNVType type = CNVType.GAINLOSS;
                    	type = cnvFeatures[3].equals("gain/loss") ? type : 
                    		CNVType.valueOf(cnvFeatures[3].toUpperCase());
                    	Region cnvEvent = new Region(chr, cnvFeatures[1], cnvFeatures[2], type);
                    	sampleCNVEvents.add(cnvEvent);
                    }
                    conradCNVEventsMap.put(caseSampleName, sampleCNVEvents);
                }
                Set<Region> predictedCNVEventsConrad = new HashSet<Region>();
                Set<Region> predictedCNVEventsTool = new HashSet<Region>();
                for (String sampleResultCNVEventLine : sampleResultCNVEventsLines) {
                	String[] cnvFeatures = sampleResultCNVEventLine.split("\\s+");
                	String chr = cnvFeatures[0].equals("1") ? "chr" + cnvFeatures[0] : 
                		cnvFeatures[0];
                	Float cnvRatio = Float.parseFloat(cnvFeatures[cnvFeatures.length - 1]);
                	CNVType type = cnvRatio > seqcnvGainCutoff ? CNVType.GAIN : CNVType.LOSS;
                	Region cnvEvent = new Region(chr, cnvFeatures[1], cnvFeatures[2], type);
                	conradCNVEventsMap.get(caseSampleName).stream().forEach(conradCNVEvent -> {
                		if (conradCNVEvent.isOverlappedWithType(cnvEvent)) {
                			predictedCNVEventsTool.add(cnvEvent);
                			predictedCNVEventsConrad.add(conradCNVEvent);
                		}
                	});
                }
                Map<String, List<String>> sampleOutput = formattedOutputResult.get(mergedControlName);
                if (sampleOutput == null) {
                	sampleOutput = new LinkedHashMap<String, List<String>>();
                }
                List<String> sampleOutputFeature = new ArrayList<String>();
                sampleOutputFeature.add(String.valueOf(sampleResultCNVEventsLines.size()));
                sampleOutputFeature.add(String.valueOf(conradCNVEventsMap.get(caseSampleName).size()));
                sampleOutputFeature.add(String.valueOf(predictedCNVEventsConrad.size()));
                sampleOutputFeature.add(String.valueOf(sampleResultCNVEventsLines.size()
                		- predictedCNVEventsTool.size()));
                sampleOutputFeature.add(String.format("%.0f%%", (float) 
                		predictedCNVEventsTool.size() * 100 / sampleResultCNVEventsLines.size()));
                sampleOutputFeature.add(String.format("%.0f%%", (float) 
                		predictedCNVEventsConrad.size() * 100 / conradCNVEventsMap.get(caseSampleName).size()));
                sampleOutput.put(caseSampleName, sampleOutputFeature);
                formattedOutputResult.put(mergedControlName, sampleOutput);
            }
        }
        for (Map.Entry<String, Map<String, List<String>>> entry 
        		: formattedOutputResult.entrySet()) {
        	System.out.println(entry.getKey());
        	List<List<String>> sampleLists = new ArrayList<List<String>>();
        	sampleLists.add(entry.getValue().get("NA10847"));
        	sampleLists.add(entry.getValue().get("NA18967"));
        	sampleLists.add(entry.getValue().get("NA18973"));
        	sampleLists.add(entry.getValue().get("NA18981"));
        	sampleLists.add(entry.getValue().get("NA19131"));
        	int i = (int) sampleLists.stream().filter(sampleList -> sampleList != null 
        			&& !sampleList.isEmpty()).count();
        	if (i != sampleLists.size()) {
        		sampleLists = sampleLists.stream().filter(sampleList -> sampleList != null 
        				&& !sampleList.isEmpty()).collect(Collectors.toList());
        	}
        	Optional<List<String>> anyNoneEmptyList = sampleLists.stream()
        			.filter(sampleList -> sampleList != null && !sampleList.isEmpty()).findAny();
        	int featureSize = 0;
        	if (anyNoneEmptyList.isPresent()) {
        		featureSize = anyNoneEmptyList.get().size();
        	}
        	for (int j = 0; j < featureSize; j++) {
        		int column = 0;
            	for (List<String> formattedList : sampleLists) {
    				System.out.print(formattedList.get(j));
    				column++;
    				if (column == sampleLists.size()) {
    					System.out.println();
    				}
    				else {
    					System.out.print("\t");
    				}
            	}
        	}
        	
        }
        
    }
    
    /**
     * Multi-case samples。使用merged control。先把case和所有merged control的sample的正常region取出来，
     * 再把case中的非正常region但是所有merged control都是正常region的region取出来。所有这些region组成了
     * 检测CNV用的bed file。先用这个bed file对SeqCNV结果过滤，所有不与该bed file中的region重叠的那些CNV结果
     * 都扔掉。然后根据那些过滤后的结果，来计算precision和recall。
     */
    public static void analyzeSeqCNVUsingDetectingBedWithMergedControl() {
    	
    }
    
    /**
     * One case sample。使用随机几个sample作merge（merge之后的bam当做control），与一个随机的sample当做case。
     * 先把case和所有merged control的sample的正常region取出来，再把case中的非正常region但是所有merged control
     * 都是正常region的region取出来。所有这些region组成了检测CNV用的bed file。先用这个bed file对SeqCNV结果过滤，
     * 所有不与该bed file中的region重叠的那些CNV结果都扔掉
     * {@link OverlapAnalysisUsingMergedControl#filterSeqCNVResultWithDetectingBed}}。
     * 然后根据那些过滤后的结果，来计算precision和recall。
     * 
     * @param filteredCaseResultFolderPath case sample result path filtered with bed.
     * @param mergedControlCNVFilePath
     * @throws IOException
     */
    public static void analyzeSeqCNVOneSampleUsingDetectingBedWithMergedControl(
    		Path filteredCaseResultFolderPath, Path mergedControlCNVFilePath) throws IOException {
    	String caseSampleName = filteredCaseResultFolderPath.getFileName().toString();
        Path caseResultFilePath = filteredCaseResultFolderPath.resolve("report")
            .resolve("filtered_CNV_report.txt");
        List<String> sampleResultCNVEventsLines = Files.readAllLines(caseResultFilePath)
        		.stream().filter(line -> {
        	String[] cnvFeatures = line.split("\\s+");
        	long cnvStartPos = Long.parseLong(cnvFeatures[1]); 
        	long cnvEndPos = Long.parseLong(cnvFeatures[2]); 
        	Float cnvRatio = Float.parseFloat(cnvFeatures[cnvFeatures.length - 1]);
        	return (cnvRatio > seqcnvGainCutoff || cnvRatio < seqcnvLossCutoff) 
        			&& (cnvEndPos - cnvStartPos) > 1000;
        }).collect(Collectors.toList());

        
        List<String> sampleCNVEventLines = Files.readAllLines(mergedControlCNVFilePath);
        sampleCNVEventLines.remove(0);  // drop header
        
        // conrad normal regions
        List<Region> normalRegions = sampleCNVEventLines.stream().filter(line -> {
        	return line.indexOf("normal") != -1;
        }).map(line -> {
        	String[] cnvFeatures = line.split("\\s+");
        	String chr = cnvFeatures[1].length() < 3 ? "chr" + cnvFeatures[1] : cnvFeatures[1];
        	// here we use GAINLOSS to stand for "normal"
        	return new Region(chr, cnvFeatures[2], cnvFeatures[3], CNVType.GAINLOSS); 
        }).collect(Collectors.toList());
        
        // conrad cnv regions
        List<Region> cnvRegions = sampleCNVEventLines.stream().filter(line -> {
        	return line.indexOf("normal") == -1;
        }).map(line -> {
        	String[] cnvFeatures = line.split("\\s+");
        	String chr = cnvFeatures[1].length() < 3 ? "chr" + cnvFeatures[1] : cnvFeatures[1];
        	Float cnvRatio = Float.parseFloat(cnvFeatures[cnvFeatures.length - 1]);
        	CNVType type = line.indexOf("gain") != -1 ? CNVType.GAIN : CNVType.LOSS;
        	return new Region(chr, cnvFeatures[2], cnvFeatures[3], type); 
        }).collect(Collectors.toList());
            
        
        Set<Region> predictedCNVEventsConrad = new HashSet<Region>();
        Set<Region> predictedCNVEventsTool = new HashSet<Region>();
        sampleResultCNVEventsLines.stream().forEach(sampleResultCNVEventLine -> {
        	String[] cnvFeatures = sampleResultCNVEventLine.split("\\s+");
        	String chr = cnvFeatures[0].equals("1") ? "chr" + cnvFeatures[0] : 
        		cnvFeatures[0];
        	Float cnvRatio = Float.parseFloat(cnvFeatures[cnvFeatures.length - 1]);
        	CNVType type = cnvRatio > seqcnvGainCutoff ? CNVType.GAIN : CNVType.LOSS;
        	Region cnvEvent = new Region(chr, cnvFeatures[1], cnvFeatures[2], type);
        	cnvRegions.stream().forEach(conradCNVEvent -> {
        		if (conradCNVEvent.isOverlappedWithType(cnvEvent)) {
        			predictedCNVEventsTool.add(cnvEvent);
        			predictedCNVEventsConrad.add(conradCNVEvent);
        		}
        	});
        });
        
        Set<Region> falsePositiveCNVEvents = new HashSet<Region>();
        sampleResultCNVEventsLines.stream().forEach(sampleResultCNVEventLine -> {
        	String[] cnvFeatures = sampleResultCNVEventLine.split("\\s+");
        	String chr = cnvFeatures[0].equals("1") ? "chr" + cnvFeatures[0] : 
        		cnvFeatures[0];
        	Float cnvRatio = Float.parseFloat(cnvFeatures[cnvFeatures.length - 1]);
        	CNVType type = cnvRatio > seqcnvGainCutoff ? CNVType.GAIN : CNVType.LOSS;
        	Region resultCNVEvent = new Region(chr, cnvFeatures[1], cnvFeatures[2], type);
        	normalRegions.stream().forEach(conradNormalCNVEvent -> {
        		if (conradNormalCNVEvent.isOverlapped(resultCNVEvent)) {
        			for (Region cnvEvent : cnvRegions) {
        				if (resultCNVEvent.isOverlappedWithType(cnvEvent)) {
        					return;
        				}
        			}
        			falsePositiveCNVEvents.add(conradNormalCNVEvent);
        		}
        	});
        });
        
        long tp = predictedCNVEventsConrad.size();
        long fn = cnvRegions.size() - tp;
        long fp = falsePositiveCNVEvents.size();
        long tn = normalRegions.size() - fp;
        
        System.out.println("control: " + mergedControlCNVFilePath.getFileName().toString());
        System.out.println("case: " + caseSampleName);
        System.out.println("seqcnv predicted cnv count: " + sampleResultCNVEventsLines.size());
        System.out.println("conrad normal region count: " + normalRegions.size());
        System.out.println("conrad cnv region count: " + cnvRegions.size());
        System.out.println("precision: " + String.format("%.0f%%", ((double) 
        		predictedCNVEventsTool.size() * 100 / sampleResultCNVEventsLines.size())));
        System.out.println("recall: " + String.format("%.0f%%", ((double) 
        		predictedCNVEventsConrad.size() * 100 / cnvRegions.size())));
        System.out.println("TruePositive: " + tp);
        System.out.println("FalseNegative: " + fn);
        System.out.println("FalsePositive: " + fp);
        System.out.println("TrueNegative: " + tn);
        System.out.println("FPR=FP/(FP+TN): " + String.format("%.0f%%", ((double) fp 
        		* 100 / (fp + tn))));
    }
    
    /**
     * filter seqcnv result with bed 
     * 
     * @param resultFilePath
     * @param filteredResultFilePath
     * @param bedFilePath
     * @return filtered seqcnv result filepath
     * @throws IOException 
     */
    public static void filterSeqCNVResultWithDetectingBed(Path resultFilePath, 
    		Path filteredResultFilePath, Path bedFilePath) throws IOException {
    	if (!Files.exists(filteredResultFilePath)) {
    		Files.createFile(filteredResultFilePath);
    	}
    	List<String> bedLines = Files.readAllLines(bedFilePath);
    	List<Region> bedRegions = bedLines.stream().map(bedLine -> {
    		String[] chrFeature = bedLine.split("\\s+");
    		String chr = chrFeature[0].length() < 3 ? "chr" + chrFeature[0] : chrFeature[0];
    		return new Region(chr, chrFeature[1], chrFeature[2]);
    	}).collect(Collectors.toList());
    	List<String> cnvResultLines = Files.readAllLines(resultFilePath);
    	List<String> filteredCnvResultLines = cnvResultLines.stream().filter(cnvLine -> {
    		String[] chrFeature = cnvLine.split("\\s+");
    		String chr = chrFeature[0].length() < 3 ? "chr" + chrFeature[0] : chrFeature[0];
    		Region cnvRegion = new Region(chr, chrFeature[1], chrFeature[2]);
    		for (Region bedRegion : bedRegions) {
    			if (cnvRegion.isOverlapped(bedRegion))
    				return true;
    		}
    		return false;
    	}).collect(Collectors.toList());
		Files.write(filteredResultFilePath, filteredCnvResultLines);
    }

}
