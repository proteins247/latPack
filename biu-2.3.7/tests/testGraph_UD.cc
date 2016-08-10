

#include <iostream>
#include <biu/Graph_UD.hh>


int main(int argc, char **argv) {

	using namespace biu;

	std::cout	<<"========== BEGIN OF TEST =================\n"
				<<" test class : Graph_UD\n"
				<<"==========================================" <<std::endl;
				
		// erzeugen	
	std::cout <<"=> creating biu::Graph_UD(5)" <<std::endl;
	Graph_UD graph = Graph_UD(5);
	
	std::cout <<"->nodes() = " <<graph.nodes() <<std::endl;
	std::cout <<"->edges() = " <<graph.edges() <<std::endl;
	
	std::vector<size_t> conID;
	std::cout <<"->connectComponents = " <<graph.connectedComponents(conID) <<std::endl;

	std::cout	<<"->addEdge(0,4)\n"	
				<<"->addEdge(0,3)\n"
				<<"->addEdge(4,3)\n"
				<<"->addEdge(2,1)" <<std::endl;
	graph.addEdge(0,4);	
	graph.addEdge(0,3);
	graph.addEdge(4,3);
	graph.addEdge(2,1);
	
	std::cout <<"->printDOT(cout) : \n";
	graph.printDOT(std::cout);
	std::cout <<std::endl;
	
	std::cout <<"->edges() = " <<graph.edges() <<std::endl;
	std::cout <<"->connectComponents = " <<graph.connectedComponents(conID) <<std::endl;
	std::cout <<"  labels = ";
	for (size_t i = 0; i< conID.size(); i++)
		std::cout <<"  "<<i <<"(" <<conID[i] <<")";
	std::cout <<std::endl;

	std::cout <<"->setNodeNumber(4)" <<std::endl;
	graph.setNodeNumber(4);	
	std::cout <<"->nodes() = " <<graph.nodes() <<std::endl;
	std::cout <<"->edges() = " <<graph.edges() <<std::endl;
	std::cout <<"->printDOT(cout) : \n";
	graph.printDOT(std::cout);
	std::cout <<std::endl;
	
	// test function remEdge
	std::cout <<"graph.remEdge(0,3)" << std::endl;
	graph.remEdge(0,3);
	std::cout <<"->printDOT(cout) : \n";
	graph.printDOT(std::cout);
	std::cout << std::endl;
	
	// test fucntion isEdge
	std::cout	<< "graph.isEdge(2,1) = " << graph.isEdge(2,1) << std::endl;
	
	// test function degree
	std::cout	<< "graph.degree(1) = " << graph.degree(1) << std::endl;

	std::cout	<<"============= END OF TEST ================\n" <<std::endl;
	
	return 0;
}
