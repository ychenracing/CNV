package wes;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
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

    public static void main(String[] args) throws IOException {
    	seqcnv();
    }

    public static void seqcnv() throws IOException {
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
                	CNVType type = cnvRatio > 1.4 ? CNVType.GAIN : CNVType.LOSS;
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
                sampleOutputFeature.add(String.format("%.0f%%", (float)predictedCNVEventsTool.size() * 100 / 
                		sampleResultCNVEventsLines.size()));
                sampleOutputFeature.add(String.format("%.0f%%", (float)predictedCNVEventsConrad.size() 
                		* 100/ conradCNVEventsMap.get(caseSampleName).size()));
                
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

}
