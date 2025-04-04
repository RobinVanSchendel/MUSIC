

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.zip.GZIPOutputStream;

import org.jcvi.jillion.core.datastore.DataStoreProviderHint;
import org.jcvi.jillion.core.util.iter.StreamingIterator;
import org.jcvi.jillion.trace.fastq.FastqDataStore;
import org.jcvi.jillion.trace.fastq.FastqFileDataStoreBuilder;
import org.jcvi.jillion.trace.fastq.FastqQualityCodec;
import org.jcvi.jillion.trace.fastq.FastqRecord;
import org.jcvi.jillion.trace.fastq.FastqWriter;
import org.jcvi.jillion.trace.fastq.FastqWriterBuilder;

public class DemultiplexMarcoLargeScreen {
	public static String unmatchedTriplet = "XXX";

	public static void main(String[] args) {
		
		//File dir = new File("Z:\\Datasets - NGS, UV_TMP, MMP\\Targeted Sequencing\\Hartwig\\LUMC-001-104-Marco_CRISPR_screen\\iSeq\\data\\");
		File dir = new File("Z:\\Datasets - NGS, UV_TMP, MMP\\Targeted Sequencing\\Hartwig\\GenomeScan105883\\RawMerged\\MB");
		File dirOut = Paths.get(dir.getAbsolutePath(),"out").toFile();
		if(!dirOut.exists()) {
			dirOut.mkdir();
		}
		String endOfBarCodeString = "ggtgtttcgtccttt".toUpperCase();
		String startOfBarCodeString = "gcatagctcttaaac".toUpperCase();
		int startOfBarCodeStringPos = 35;
		int endOfBarCodeStringPos = 54;
		boolean write = true;
		boolean countNonMatched = false;
		boolean addBarcodeToFQ = true;
		boolean separateInBins = true;
		
		boolean outputR2 = false;
		
		ArrayList<PairedEnd> files = getPairedEnd(dir);
		//HashMap<String, String> barcodes = createHashMap("104596_barcodes_MB.txt");
		String yusa = "yusa.txt";
		HashMap<String, String> barcodes = createHashMap(yusa);
		HashMap<String, String> sgRNAToBins = createHashMapSGRNAs(yusa);
		try {
			for(PairedEnd pe: files) {
				if(!pe.getR1().getName().contains("non-t")) {
					//continue;
				}
				HashMap<String, FastqWriter> hmWriterR1 = new HashMap<String, FastqWriter>();
				HashMap<String, FastqWriter> hmWriterR2 = new HashMap<String, FastqWriter>();
				ArrayList<File> fileList = new ArrayList<File>();
				
				FastqDataStore datastoreR1 = new FastqFileDataStoreBuilder(pe.getR1())
						.qualityCodec(FastqQualityCodec.SANGER)
	                    .hint(DataStoreProviderHint.ITERATION_ONLY)
	                    .build();
				FastqDataStore datastoreR2 = new FastqFileDataStoreBuilder(pe.getR2())
						.qualityCodec(FastqQualityCodec.SANGER)
	                    .hint(DataStoreProviderHint.ITERATION_ONLY)
	                    .build();
				
				StreamingIterator<FastqRecord> iterR1 = datastoreR1.iterator();
				StreamingIterator<FastqRecord> iterR2 = datastoreR2.iterator();
				
				int nr = 0;
				int hit = 0;
				int nonhit = 0;
				int notFoundSG = 0;
				int hitPosition = 0;
				
				HashMap<String, Integer> barcodeHits = new HashMap<String,Integer>();
				while(iterR1.hasNext()) {
					FastqRecord fqR1 = iterR1.next();
					FastqRecord fqR2 = iterR2.next();
					//String R1S = fqR1.getNucleotideSequence().toString();
					String R2S = fqR2.getNucleotideSequence().toString();
					
					//System.out.println(R1S.substring(0, endOfBarCode));
					int start = R2S.indexOf(startOfBarCodeString);
					int end = R2S.indexOf(endOfBarCodeString);
					if(start>=0 && end>=0 && end > (start+startOfBarCodeString.length())) {
						start += startOfBarCodeString.length();
						String barcode = R2S.substring(start, end-1);
						
						if(barcodes.containsKey(barcode)) {
							hit++;
						}
						//try position
						else {
							barcode = R2S.substring(startOfBarCodeStringPos,endOfBarCodeStringPos);
							if(barcodes.containsKey(barcode) ) {
								hitPosition++;
							}
						}
						
						if(barcodes.containsKey(barcode) ) {
							String foundBarcode = barcodes.get(barcode);
							//overwrite barcode if required
							String barcodeLookup = foundBarcode;
							if(addBarcodeToFQ && separateInBins) {
								foundBarcode = sgRNAToBins.get(foundBarcode);
							}
							else if(addBarcodeToFQ) {
								foundBarcode = "therecanonlybeone";
							}
							if(!barcodeHits.containsKey(foundBarcode)){
								barcodeHits.put(foundBarcode, 0);
							}
							barcodeHits.put(foundBarcode, barcodeHits.get(foundBarcode)+1);
							
							//Write to file
							if(!hmWriterR1.containsKey(foundBarcode)) {
								String fileName = pe.getR1().getName().replace("_R1_001.fastq.gz", "");
								fileName = fileName.replace("R1.fastq.gz", "");
								if(write) {
									File R1 = Paths.get(dirOut.getAbsolutePath(),fileName+"_"+foundBarcode+"_R1.fastq").toFile(); 
									hmWriterR1.put(foundBarcode, new FastqWriterBuilder(R1).build());
									fileList.add(R1);
									if(outputR2) {
										File R2 = Paths.get(dirOut.getAbsolutePath(),fileName+"_"+foundBarcode+"_R2.fastq").toFile();
										hmWriterR2.put(foundBarcode, new FastqWriterBuilder(R2).build());
										fileList.add(R2);
									}
									
									
									
								}
								else {
									hmWriterR1.put(foundBarcode,null);
									hmWriterR2.put(foundBarcode,null);
								}
							}
							if(write) {
								if(addBarcodeToFQ) {
									FastqRecord fqTestR1 = fqR1.toBuilder().comment("BC:"+barcodeLookup).build();
									hmWriterR1.get(foundBarcode).write(fqTestR1);
									if(outputR2) {
										FastqRecord fqTestR2 = fqR2.toBuilder().comment("BC:"+barcodeLookup).build();
										hmWriterR2.get(foundBarcode).write(fqTestR2);
									}
									
								}
								else {
									hmWriterR1.get(foundBarcode).write(fqR1);
									hmWriterR2.get(foundBarcode).write(fqR2);
								}
							}
						}
						else {
							if(countNonMatched) {
								if(!barcodeHits.containsKey(barcode)){
									barcodeHits.put(barcode, 0);
								}
								barcodeHits.put(barcode, barcodeHits.get(barcode)+1);
							}
							//System.out.println("nonhit\t"+pe.getR1().getName()+"\t"+barcodeAlt+"\t"+Utils.reverseComplement(barcodeAlt));
							nonhit++;
						}
						
					}
					else {
						//System.out.println(start+"\t"+end+"\t"+R2S);
						notFoundSG++;
					}
					
					nr++;
					if(nr%100000==0) {
						System.out.println("Already processed "+nr+" records");
					}
				}
				System.out.println(pe.getR1().getName()+"\thit "+hit+" hit position "+hitPosition+" non-hit: "+nonhit+" not-sgRNAsurround: "+notFoundSG);
				iterR1.close();
				iterR2.close();
				if(write) {
					for(String key: hmWriterR1.keySet()) {
						hmWriterR1.get(key).close();
					}
					for(String key: hmWriterR2.keySet()) {
						hmWriterR2.get(key).close();
					}
					//compress and delete
					for(File f: fileList) {
						System.out.println("Compressing "+f.getName());
						compressGZIP(f);
						f.delete();
					}
				}
				
			}
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (RuntimeException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static HashMap<String, String> createHashMapSGRNAs(String string) {
		HashMap<String, String> hm = new HashMap<String, String>();
		int binNameColumn = 7;
		int idColumn = 2;
		try {
			Scanner s = new Scanner(new File(string));
			while(s.hasNextLine()) {
				String line = s.nextLine();
				String[] parts = line.split("\t");
				String name = parts[idColumn];
				String bin = parts[binNameColumn];
				hm.put(name, bin);
			}
			s.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(0);
		}
		System.out.println("Added "+hm.size()+" barcodes");
		return hm;
	}
	private static Map<String, String> sortByComparator(Map<String, String> unsortMap, final boolean order)
    {

        List<Entry<String, String>> list = new LinkedList<Entry<String, String>>(unsortMap.entrySet());

        // Sorting the list based on values
        Collections.sort(list, new Comparator<Entry<String, String>>()
        {
            public int compare(Entry<String, String> o1,
                    Entry<String, String> o2)
            {
                if (order)
                {
                    return o1.getValue().compareTo(o2.getValue());
                }
                else
                {
                    return o2.getValue().compareTo(o1.getValue());

                }
            }
        });

        // Maintaining insertion order with the help of LinkedList
        Map<String, String> sortedMap = new LinkedHashMap<String, String>();
        for (Entry<String, String> entry : list)
        {
            sortedMap.put(entry.getKey(), entry.getValue());
        }

        return sortedMap;
    }

	private static ArrayList<PairedEnd> getPairedEnd(File f) {
		ArrayList<PairedEnd> al = new ArrayList<PairedEnd>();
		if(f.isDirectory()) {
			for(File R1: f.listFiles()) {
				if(R1.getName().endsWith("R1_001.fastq.gz")) {
					String str = R1.getAbsolutePath().replace("R1_001.fastq.gz", "R2_001.fastq.gz");
					File R2 = new File(str);
					if(R2.exists()) {
						PairedEnd pe = new PairedEnd(R1, R2);
						al.add(pe);
					}
				}
				else if(R1.getName().endsWith("_1.fastq")) {
					String str = R1.getAbsolutePath().replace("_1.fastq", "_2.fastq");
					File R2 = new File(str);
					if(R2.exists()) {
						PairedEnd pe = new PairedEnd(R1, R2);
						al.add(pe);
					}
				}
				else if(R1.getName().endsWith("R1.fastq.gz")){
					String str = R1.getAbsolutePath().replace("R1.fastq.gz", "R2.fastq.gz");
					File R2 = new File(str);
					String str3 = R1.getAbsolutePath().replace("R1.fastq.gz", "R3.fastq.gz");
					File R3 = new File(str3);
					if(R3.exists()) {
						PairedEnd pe = new PairedEnd(R1, R3);
						al.add(pe);
					}
					else if(R2.exists()) {
						PairedEnd pe = new PairedEnd(R1, R2);
						al.add(pe);
					}
				}
			}
		}
		return al;
		
	}

	private static HashMap<String, String> createHashMap(String string) {
		HashMap<String, String> hm = new HashMap<String, String>();
		System.out.println("Taking column 2 and 1 for barcode and ID");
		int idColumn = 2;
		int sgColumn = 3;
		boolean takeRevCom = true;
		boolean removeAAAC = false;
		try {
			Scanner s = new Scanner(new File(string));
			while(s.hasNextLine()) {
				String line = s.nextLine();
				String[] parts = line.split("\t");
				String revCom = parts[sgColumn].toUpperCase();
				if(takeRevCom) {
					revCom = Utils.reverseComplement(revCom);
				}
				//FastqWriter fw = new FastqWriterBuilder(outFile).build();
				String id = parts[idColumn].replaceAll("[^a-zA-Z0-9\\.\\-]", "_");
				
				if(removeAAAC) {
					revCom = revCom.substring(4);
				}
				//System.out.println("Adding ["+revCom+"] "+revCom.length());
				hm.put(revCom, id);
			}
			s.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(0);
		}
		System.out.println("Added "+hm.size()+" barcodes");
		return hm;
	}

	private static ArrayList<File> getGZFiles(File f) {
		ArrayList<File> al = new ArrayList<File>();
		if(f.isDirectory()) {
			for(File file: f.listFiles()) {
				if(file.getName().endsWith(".fastq.gz")) {
					al.add(file);
				}
			}
		}
		return al;
	}
	public static void compressGZIP(File input) throws IOException {
		File output = new File(input.getAbsolutePath()+".gz");
        try (GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(output))){
            try (FileInputStream in = new FileInputStream(input)){
                byte[] buffer = new byte[1024];
                int len;
                while((len=in.read(buffer)) != -1){
                    out.write(buffer, 0, len);
                }
            }
        }
    }

}
