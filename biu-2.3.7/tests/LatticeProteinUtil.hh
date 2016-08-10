#ifndef LATTICEPROTEINUTIL_HH_
#define LATTICEPROTEINUTIL_HH_

#include "biu/assertbiu.hh"
#include "biu/LatticeModel.hh"
#include "biu/LatticeDescriptorSQR.hh"

#include <sstream>

	  /*! Converts an absolute move sequence into an ASCII representation
	   * of its 2D coordinates.
	   * 
	   * @param moves the move sequence to convert
	   * 
	   * @return the 2D coordinates string representation
	   * 
	   */
	std::string
	toASCII( const biu::MoveSequence & moves )
	{
		assertbiu(moves.size() > 0, "no moves given");
		
		using namespace biu;
		
		LatticeDescriptorSQR ld;
		LatticeModel lattice(&ld);
		
		  // convert to points
		IPointVec p = lattice.absMovesToPoints( moves );

		  // get min/max dimension
		int xMin = 0;
		int xMax = 0;
		int yMin = 0;
		int yMax = 0;
		
		for (size_t i=0; i<p.size();i++) {
			xMin = std::min(p[i].getX(),xMin);
			xMax = std::max(p[i].getX(),xMax);
			yMin = std::min(p[i].getY(),yMin);
			yMax = std::max(p[i].getY(),yMax);
			
//			std::cout <<" (" <<p[i] <<")\n";
		}
		
		  // create output panel
		std::vector < std::string > out( 2*(xMax-xMin+1), std::string( 2*(yMax-yMin+1), ' ') );
//		std::vector < std::string > out( 2*(yMax-yMin+1), std::string( 2*(xMax-xMin+1), ' ') );
		
		  // fill output panel
		for (size_t i=0; i<p.size(); i++) {
			  // calculate relative coordinates in output
			int xPos = 2 * (xMax-p[i].getX());
			int yPos = 2 * (p[i].getY()-yMin);
			  // print monomer
			out[xPos][yPos] = 'O';
//			out[yPos][xPos] = 'O';
			
			  // set backbone
			if (i>0) {
				int xDiff = p[i].getX() - p[i-1].getX();
				int yDiff = p[i].getY() - p[i-1].getY();
				
				if (yDiff != 0) {
					out[xPos][yPos-yDiff] = '-';
				} else {
					out[xPos+xDiff][yPos] = '|';
				}
//				if (xDiff != 0) {
//					out[yPos][xPos-xDiff] = '-';
//				} else {
//					out[yPos-yDiff][xPos] = '|';
//				}
			}
		}
		
		  // create return string
		std::stringstream ret;
		for (size_t i=0; (i+1)<out.size(); i++) {
			ret <<" "<< out[i] <<"\n";
		}
		
		return ret.str();
	}

#endif /*LATTICEPROTEINUTIL_HH_*/
