
public class Utils {
	public static String reverseComplement(String dna){
		StringBuffer dnaRev = new StringBuffer(dna).reverse();
		String revCom = "";
		for(char c: dnaRev.toString().toCharArray()){
			switch(c){
				case 'a':
					revCom += 't';
					break;
				case 'A':
					revCom += 'T';
					break;
				case 't':
					revCom += 'a';
					break;
				case 'T':
					revCom += 'A';
					break;
				case 'c':
					revCom += 'g';
					break;
				case 'C':
					revCom += 'G';
					break;
				case 'g':
					revCom += 'c';
					break;
				case 'G':
					revCom += 'C';
					break;
				case 'N':
					revCom += 'N';
					break;
				case 'Y':
					revCom += 'Y';
					break;
				default:
					System.err.println("Can't complement "+c);
					System.err.println("Can't complement "+dna);
			}
		}
		return revCom;
	}
}
