

#include <iostream>
#include <biu/Alphabet.hh>


int main(int argc, char **argv) {

	using namespace biu;

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test class : Alphabet\n"
				<<"==========================================" <<std::endl;
				
		// erzeugen	
	std::cout <<"=> creating biu::Alphabet(HPNX,1)" <<std::endl;
	Alphabet alph = Alphabet("HPNX",1);
	
		// getAlphabetSize()
	std::cout <<"->getAlphabetSize() = " <<alph.getAlphabetSize() <<std::endl;
	
		// getElementLength()
	std::cout	<<"->getElementLength() = " << alph.getElementLength() 
				<<std::endl;
	
		// isAlphabetString() == false
		std::string str = "HPNXY";
	std::cout	<<"->isAlphabetString("<<str<<")  = " 
				<< (alph.isAlphabetString(str)?"true":"false")
				<<std::endl;

		// isAlphabetString() == true
	str = "HPPPHPHPHNNNXXXPPHHXXX";
	std::cout	<<"-> " <<str <<" = str" <<std::endl;
	std::cout	<<"->isAlphabetString(str)  = " 
				<< (alph.isAlphabetString(str)?"true":"false")
				<<std::endl;

		// getElement()	
	std::cout	<<"->getElement(H)  = " << alph.getElement("H")
				<<std::endl;

		// getString( getElement( X ) ) == X
	std::cout	<<"->getString(getElement(H))  = " 
				<< alph.getString(alph.getElement("H"))
				<<std::endl;

		// getString( getSequence( X ) ) == X
	std::cout	<<"-> (getString(getSequence(str) == str)  = " 
				<< (alph.getString(alph.getSequence(str))==str?"true":"false")
				<<std::endl;
				
		// getIndex()
	std::cout	<<"->getIndex(X)  = " << alph.getIndex("X")
				<<std::endl;

		// getIndex(ELEM)
	std::cout	<<"->getIndex(getElement(X))  = " 
				<< alph.getIndex(alph.getElement("X"))
				<<std::endl;

		// getIndex()
	std::cout	<<"->getString(getElement(getIndex(X)))  = " 
				<< alph.getString(alph.getElement(alph.getIndex("X")))
				<<std::endl;

		// compress / decompress
	std::cout	<<"-> (getString(decompress(compress(getSequence(str)))) == str) = "
				<< ( alph.getString(alph.decompress(alph.compress(alph.getSequence(str)), str.size())) == str ? "true" : "false")
				<<std::endl;
	std::cout	<<"-> (getString(decompress(compressS(str))) == str) = "
				<< ( alph.getString(alph.decompress(alph.compressS(str), str.size())) == str ? "true" : "false")
				<<std::endl;
	std::cout	<<"-> (decompressS(compressS(str)) == str) = "
				<< ( alph.decompressS(alph.compressS(str), str.size()) == str ? "true" : "false")
				<<std::endl;

	
/*
		// 
	std::cout	<<"->  = "
				<<lp.
				<<std::endl;

*/

	std::cout	<<"============= END OF TEST ================\n" <<std::endl;
	
	return 0;
}
