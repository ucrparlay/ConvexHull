#include "get_time.h"
#include "ConvexHull.h"

int main(int argc, char** argv) {
	
	if (argc != 3) {
		cout << "Usage: ./qsort [num_elements] [type of input]" << endl;
		return 0;
	}
	int n = atoi(argv[1]);
	int type_of_input = atoi(argv[2]);

	point* p = new point[n];
	
	if(type_of_input == 0){
		RandmPoint(p, n);
	}
	else if (type_of_input == 1) {
		RandmOnCircle(p, n);
	}

	init(p,n);

	timer timer1;
	facet facet1, facet2;
	timer1.start();
	if (start(p[0], p[1], p[2])) {
		ProcessRidge(facet1 = { p[1],p[0] }, p[0], facet2 = { p[0],p[2] });
		ProcessRidge(facet1 = { p[2],p[1] }, p[1], facet2 = { p[1],p[0] });
		ProcessRidge(facet1 = { p[0],p[2] }, p[2], facet2 = { p[2],p[1] });
	}
	else {
		ProcessRidge(facet1 = { p[2],p[0] }, p[0], facet2 = { p[0],p[1] });
		ProcessRidge(facet1 = { p[0],p[1] }, p[1], facet2 = { p[1],p[2] });
		ProcessRidge(facet1 = { p[1],p[2] }, p[2], facet2 = { p[2],p[0] });
	}
	timer1.stop();

	cout << "total time: " << timer1.get_total() << endl;
	cout << "output size: " << H.size() << endl;


	return 0;
}
