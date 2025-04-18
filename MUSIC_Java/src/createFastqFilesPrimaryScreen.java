

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.Scanner;
import java.util.zip.GZIPOutputStream;

public class createFastqFilesPrimaryScreen {
	
	public static final String LEFTSEQ = "AACTTGCTATGCTGTTTCCAGCATAGCTCTTAAAC";
	public static final String RIGHTSEQ = "CGGTGTTTCGTCCTTTCCACAAGATATATAAAGCCAAGAAATCGAAATACTTTCAAGTTACGGTAAGCATATGATAGTCCATTTTAAAACATAATTT";

	public static void main(String[] args) {
		File f = new File(args[0]);
		File lib = new File(args[1]);
		lib = new File("yusa_orig.txt");
		HashMap<String, String> barcodeToSeq = null;
		
		try {
			barcodeToSeq = createHashMap(lib);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		
		//f = new File("Z:\\Datasets - NGS, UV_TMP, MMP\\Targeted Sequencing\\Hartwig\\LUMC-001-104-Marco_CRISPR_screen\\Raw\\Processed_Again\\MB01.txt");
		int lastDot = f.getAbsolutePath().lastIndexOf(".");
		String output = f.getAbsolutePath().substring(0, lastDot)+"_R1.fastq.gz";
		String outputR2 = f.getAbsolutePath().substring(0, lastDot)+"_R2.fastq.gz";
		File outputFR1 = new File(output);
		File outputFR2 = new File(outputR2);
		System.out.println("Input file "+f.getAbsolutePath());
		System.out.println("Output file R1\t"+outputFR1.getAbsolutePath());
		System.out.println("Output file R2\t"+outputFR2.getAbsolutePath());
		createFastQFile(f, outputFR1, outputFR2, barcodeToSeq);
	}

	private static HashMap<String, String> createHashMap(File lib) throws FileNotFoundException {
		HashMap<String, String> hm = new HashMap<String, String>();
		Scanner s = new Scanner(lib);
		int barcodeIndex = 2;
		int seqIndex = 3;
		boolean takeRC = true;
		//remove header
		if(s.hasNextLine()) {
			s.nextLine();
		}
		
		while(s.hasNextLine()) {
			String line = s.nextLine();
			String[] parts = line.split("\t");
			String seq = parts[seqIndex];
			if(takeRC) {
				seq = Utils.reverseComplement(seq);
			}
			hm.put(parts[barcodeIndex], seq);
		}
		s.close();
		return hm;
	}

	private static void createFastQFile(File in, File R1, File R2, HashMap<String, String> barcodeToSeq) {
		try {
			Scanner s = new Scanner(in);
			

			FileOutputStream fos = new FileOutputStream(R1);
			GZIPOutputStream gos = new GZIPOutputStream(fos);
			OutputStreamWriter osw = new OutputStreamWriter(gos);
			
			FileOutputStream fos2 = new FileOutputStream(R2);
			GZIPOutputStream gos2 = new GZIPOutputStream(fos2);
			OutputStreamWriter osw2 = new OutputStreamWriter(gos2);
			

			
			BufferedWriter R1w = new BufferedWriter(osw);
			BufferedWriter R2w = new BufferedWriter(osw2);
			String header = s.nextLine();
			String[] headers = header.split("\t");
			int rawC = -1;
			int idC = -1;
			int countsC = -1;
			int barcodeC = -1;
			//find the correct columns
			for(int i = 0;i< headers.length;i++) {
				if(headers[i].contentEquals("Raw")) {
					rawC = i;
				}
				else if(headers[i].contentEquals("Name")) {
					idC = i;
				}
				else if(headers[i].contentEquals("countEvents")) {
					countsC = i;
				}
				else if(headers[i].contentEquals("Barcode")) {
					barcodeC = i;
				}
				//System.out.println(i+"\t"+headers[i]);
			}
			HashMap<Integer, String> bq = new HashMap<Integer, String>();
			int lines = 0;
			while(s.hasNextLine()) {
				String line = s.nextLine();
				String[] parts = line.split("\t");
				String seq = parts[rawC];
				String ID = parts[idC];
				int counts = Integer.parseInt(parts[countsC]);
				if(seq.contains("X") ){
						seq = seq.replaceAll("X", "N");
						//System.out.println("replaced X by N");
				}
				for(int i=0;i<counts;i++) {
					String tempId = ID+":"+(i+1);
					String bqStr = getBaseQuality(seq, bq);
					String barcodeSeq = getBarcodeSeq(parts[barcodeC], barcodeToSeq);
					writeFastq(R1w,seq, tempId, bqStr);
					writeFastq(R2w,barcodeSeq, tempId, bqStr);
				}
				lines++;
				if(lines % 10000 == 0) {
					System.out.println("Already processed "+lines+" lines");
				}
			}
			s.close();
			R1w.close();
			R2w.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
	private static String getBarcodeSeq(String barcode, HashMap<String, String> barcodeToSeq) {
		String s = LEFTSEQ + barcodeToSeq.get(barcode) + RIGHTSEQ;
		return s;
	}

	//the fastq is written here based on the sequence and ID
	//basequality is set to maximum
	private static void writeFastq(BufferedWriter r1w, String seq, String id, String qual) {
		
		try {
			r1w.write("@"+id+"\n");
			r1w.write(seq+"\n");
			r1w.write("+\n");
			r1w.write(qual+"\n");
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
	}
	
	//to speed up the base quality stuff I have added a hash as seq lengths are anywhere between 20-150
	private static String getBaseQuality(String seq, HashMap<Integer, String> bqHash) {
		String bq = bqHash.get(seq.length());
		if(bq != null) {
			return bq;
		}
		StringBuffer bqB = new StringBuffer(seq.length());
		for(int i=0;i<seq.length();i++) {
			bqB.append("F");
		}
		String ret = bqB.toString();
		bqHash.put(seq.length(), ret);
		return ret;
		
	}

}
