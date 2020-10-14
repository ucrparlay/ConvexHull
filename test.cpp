#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <list>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <atomic>
#include <mutex>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "get_time.h"
#include "utilities.h"
#include "../parlaylib/include/parlay/sequence.h"
#include "../parlaylib/include/parlay/hash_table.h"


const  double PI = 3.1415926;
const  double EPSINON = 0.0000001;
int count1 = 0;
int count2 = 0;
int count3 = 0;
int count4 = 0;
int count5 = 0;
timer testtime1, testtime2, testtime3, testtime4, testtime5;


inline uint32_t hash1(uint32_t a) {
	a = (a + 0x7ed55d16) + (a << 12);
	a = (a ^ 0xc761c23c) ^ (a >> 19);
	a = (a + 0x165667b1) + (a << 5);
	a = (a + 0xd3a2646c) ^ (a << 9);
	a = (a + 0xfd7046c5) + (a << 3);
	a = (a ^ 0xb55a4f09) ^ (a >> 16);
	if (a < 0) a = -a;
	return a;
}

inline uint32_t hash2(uint32_t a){
	a = (a + 0x1a976ce5) + (a << 15);
	a = (a ^ 0x98bf7c21) ^ (a >> 13);
	a = (a + 0xc092738a) + (a << 7);
	a = (a + 0xbbc2976e) ^ (a << 9);
	a = (a + 0xf6efd123) + (a << 11);
	a = (a ^ 0x987abc01) ^ (a >> 2);
	if (a < 0) a = -a;
	return a;
}

struct point
{
	double x, y;
	int pivot,n;
	bool operator < (const point &o) const
	{
		return pivot < o.pivot;
	}
	bool operator == (const point &o) const
	{
		return pivot == o.pivot;
	}/**/
};

typedef struct facet
{
	point v1, v2;
	bool operator < (const facet &o) const
	{
		if (v1.pivot == o.v1.pivot)
			return v2.pivot < o.v2.pivot;
		return v1.pivot < o.v1.pivot;
	}
	bool operator == (const facet &o) const
	{
		return v1.pivot == o.v1.pivot && v2.pivot == o.v2.pivot;
	}
}Facet;

struct mapCkey
{
	int p1, p2;
};

struct mapCvalue
{
	parlay::sequence<point> it;
	int size;
};

int hashfacet(facet o) {
	
	return ((hash1(o.v1.pivot) ^ hash1(o.v2.pivot) << 1) % (o.v1.n*10));
}

using namespace std;

//unordered_map<facet, parlay::sequence<point>> mapC;
//unordered_map<point, facet> M;
facet* H;
point* p;
pair <int, int>* R;
pair<facet,mapCvalue>* hashC;
parlay::hashtable<parlay::hash_numeric<int>> table(10000000, parlay::hash_numeric<int>{});

int checkc = 0;
parlay::sequence<point> checkit;

bool sortfunction(point i, point j) { return (j < i); }

bool mapCinsert(facet t,mapCvalue val) {
	int i = hashfacet(t);
	while (!table.insert(i)) {
		if (hashC[i].first == t) {
			return false;
		}
		i = (i + 1) % (t.v1.n * 10);
	}
	hashC[i] = pair<facet, mapCvalue>(t, val);
	return true;
}

mapCvalue mapCfind(facet t) {
	int i = hashfacet(t);

	while (!(hashC[i].first == t)) {
		i = (i + 1) % (t.v1.n * 10);
	}
	mapCvalue out = hashC[i].second;
	
	return out;
}

bool InsertandSetR(facet t) {
	int i = hash1(t.v1.pivot) % (t.v1.n * 10);
	while (!pbbs::atomic_compare_and_swap(&R[i], pair <int, int> (0,0), pair<int, int>(t.v1.pivot, t.v2.pivot))) {
		if (R[i].first == t.v1.pivot) {
			return false;
		}
		i = (i + 1) % (t.v1.n * 10);
	}
	return true;
}

bool InsertandSetL(facet t) {
	int i = hash1(t.v2.pivot) % (t.v2.n * 10);
	while (!pbbs::atomic_compare_and_swap(&R[i], pair <int, int>(0, 0), pair<int, int>(t.v2.pivot, t.v1.pivot))) {
		if (R[i].first == t.v2.pivot) {
			return false;
		}
		i = (i + 1) % (t.v2.n * 10);
	}
	return true;
}

facet getvalueL(facet t){
	int i = hash1(t.v2.pivot) % (t.v2.n * 10);
	while (!(R[i].first == t.v2.pivot)){
		i = (i + 1) % (t.v2.n * 10);
	}
	facet out = { p[R[i].first], p[R[i].second] };
	return out;
}

facet getvalueR(facet t) {
	int i = hash1(t.v1.pivot) % (t.v1.n * 10);
	while (!(R[i].first == t.v1.pivot)) {
		i = (i + 1) % (t.v1.n * 10);
	}
	facet out = { p[R[i].second], p[R[i].first] };
	return out;
}

void RandmOnCircle(point *p, int n) {
	for (int i = 0; i < n; i++) {
		p[i].x = (double)n + n*cos((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
		p[i].y = (double)n + n*sin((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
		p[i].pivot = i;
		p[i].n = n;
	}
	return;
}

bool visible(point v, facet t) {
	if (EPSINON < (t.v1.x - v.x) * (t.v2.y - v.y) - (t.v1.y - v.y) * (t.v2.x - v.x)) return true;
	else return false;
}

bool start(point v1, point v2, point v3) {
	if (EPSINON < (v1.x - v3.x) * (v2.y - v3.y) - (v1.y - v3.y) * (v2.x - v3.x)) return true;
	else return false;
}

bool collinear(point v1, point v2, point v3) {
	if ((v1.x - v3.x) * (v2.y - v3.y) - (v1.y - v3.y) * (v2.x - v3.x) < EPSINON  && (v1.x - v3.x) * (v2.y - v3.y) - (v1.y - v3.y) * (v2.x - v3.x) > -EPSINON) return true;
	else return false;
}



void ProcessRidge(facet t1, point r, facet t2) {
	
	testtime1.start();

	parlay::sequence<point>::const_iterator t1iter;
	parlay::sequence<point>::const_iterator t2iter;
	
	auto vp = parlay::sequence<point>();
	int t1size, t2size;
	
	mapCvalue t1val = mapCfind(t1);
	
	t1iter = t1val.it.begin();
	mapCvalue t2val = mapCfind(t2);
	t2iter = t2val.it.begin();
	
	testtime1.stop();
	if ((t1val.size == 0) && (t2val.size == 0)) {
		H[t1.v1.pivot] = t1;
		count1++;
		
	}
	else if (!(t2val.size == 0) &&  !(t1val.size == 0) && *(t1val.it.begin()) == *(t2val.it.begin())) {
		count2++;
		testtime1.start();
		testtime1.stop();
	}
	else if ((t1val.size == 0) || (t2val.size != 0 && (*(t2val.it.begin()) < *(t1val.it.begin())))) {//左
		testtime2.start();
		count3++;
		parlay::sequence<point> c;
		facet t;
		t.v1 = r;
		t.v2 = *(t2val.it.begin());
		if (t1val.size == 0) {
			for (int i = 0; i < t2val.size; i++) {
				count5++;
				if (visible(*t2iter, t)) {
					vp.push_back((*t2iter));
				}
				t2iter++;
			}
		}
		else {
			int pointer1 = 0;
			int pointer2 = 0;
			while (pointer1 < t1val.size && pointer2 < t2val.size)
			{
				if (*(t1iter + pointer1) < *(t2iter + pointer2)) {
					if (visible(*(t1iter + pointer1), t)) {
						vp.push_back(*(t1iter + pointer1));
					}
					pointer1++;
				}
				else if (*(t2iter + pointer2) < *(t1iter + pointer1)) {
					if (visible(*(t2iter + pointer2), t)) {
						vp.push_back(*(t2iter + pointer2));
					}
					pointer2++;
				}
				else {
					if (visible(*(t1iter + pointer1), t)) {
						vp.push_back(*(t1iter + pointer1));
					}
					pointer1++;
					pointer2++;
				}
			}
			if (pointer1 != t1val.size) {
				for (int i = pointer1; i < t1val.size; i++) {
					if (visible(*(t1iter + i), t)) {
						vp.push_back(*(t1iter + i));
					}
				}
			}
			else if (pointer2 != t2val.size) {
				for (int i = pointer2; i < t2val.size; i++) {
					if (visible(*(t2iter + i), t)) {
						vp.push_back(*(t2iter + i));
					}
				}
			}
		}
		
		testtime2.stop();
		testtime3.start();
		mapCvalue tempval;
		tempval.it = vp;
		tempval.size = vp.size();
		mapCinsert(t, tempval);

		
		testtime3.stop();
		cilk_for (int i = 0; i < 2; i++) {
			if (i) {
				ProcessRidge(t1, r, t);
			}
			else if(!InsertandSetL(t)){
				facet oldt = getvalueL(t);
				ProcessRidge(t, t.v2, oldt);
			}
		}
	}
	else {//右
		testtime5.start();
		count4++;
		parlay::sequence<point> c;
		facet t;
		t.v1 = *(t1val.it.begin());
		t.v2 = r;
		if (t2val.size == 0) {
			for (int i = 0; i < t1val.size; i++) {
				count5++;
				if (visible(*(t1iter + i), t)) {
					vp.push_back(*(t1iter + i));
				}
			}
		}
		else {
			int pointer1 = 0;
			int pointer2 = 0;
			while (pointer1 < t1val.size && pointer2 < t2val.size)
			{
				if (*(t1iter + pointer1) < *(t2iter + pointer2)) {
					if (visible(*(t1iter + pointer1), t)) {
						vp.push_back(*(t1iter + pointer1));
					}
					pointer1++;
				}
				else if (*(t2iter + pointer2) < *(t1iter + pointer1)) {
					if (visible(*(t2iter + pointer2), t)) {
						vp.push_back(*(t2iter + pointer2));
					}
					pointer2++;
				}
				else {
					if (visible(*(t1iter + pointer1), t)) {
						vp.push_back(*(t1iter + pointer1));
					}
					pointer1++; 
					pointer2++;
				}
			}
			if (pointer1 != t1val.size) {
				for (int i = pointer1; i < t1val.size; i++) {
					if (visible(*(t1iter + i), t)) {
						vp.push_back(*(t1iter + i));
					}
				}
			}
			else if (pointer2 != t2val.size) {
				for (int i = pointer2; i < t2val.size; i++) {
					if (visible(*(t2iter + i), t)) {
						vp.push_back(*(t2iter + i));
					}
				}
			}
	}	

		mapCvalue tempval;
                tempval.it = vp;
                tempval.size = vp.size();
                mapCinsert(t, tempval);

		
		point tempr;
		facet tempt;
		testtime5.stop();
		cilk_for (int i = 0; i < 2; i++) {
			if (i) {
				ProcessRidge(t, r, t2);
			}
			else if (!InsertandSetR(t)) {
				facet oldt = getvalueR(t);
				ProcessRidge(oldt, t.v1,t );
			}
		}
	}
	return;
}

int main(int argc, char** argv) {
	
	if (argc != 3) {
		cout << "Usage: ./qsort [num_elements] [type of input]" << endl;
		return 0;
	}
	int n = atoi(argv[1]);
	int type_of_input = atoi(argv[2]);
	p = new point[n];
	R = new pair<int, int>[10 * n];
	H = new facet[n];
	hashC = new pair<facet, mapCvalue>[10 * n];
	parlay::sequence<point> hashtest;
	parlay::sequence<point>::const_iterator it;
	if(type_of_input == 0){
		for (int i = 0; i < n; i++) {
			p[i].x = (hash1(i)) % (n * 2);
			p[i].y = (hash2(i)) % (n * 2);
			p[i].pivot = i;
			p[i].n = n;
		}
	}
	else if (type_of_input == 1) {
		RandmOnCircle(p, n);
	}


	while (collinear(p[0], p[1], p[2])) {//排除初始三点共线
		int x = rand() % n;
		point temp = p[2];
		p[2] = p[x];
		p[x] = temp;
	}
	parlay::sequence<point> visiblep;
	point temp;
	facet t, facet11, facet12, facet21, facet22, facet31, facet32;
	mapCvalue tempval;
	//确定每条线中的点顺时针排列
	//
	cout << "初始C,H" << endl;
	if (start(p[0], p[1], p[2])) {
		parlay::sequence<point> visiblep1;
		t = { p[0],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep1.push_back(p[j]);
			}
		}
		sort(visiblep1.begin(), visiblep1.end());
		tempval.it = visiblep1;
                tempval.size = visiblep1.size();
		mapCinsert(t, tempval);

		parlay::sequence<point> visiblep2;
		t = { p[2],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep2.push_back(p[j]);
			}
		}
		sort(visiblep2.begin(), visiblep2.end());
		tempval.it = visiblep2;
                tempval.size = visiblep2.size();
                mapCinsert(t, tempval);

		parlay::sequence<point> visiblep3;
		t = { p[1],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep3.push_back(p[j]);
			}
		}
		sort(visiblep3.begin(), visiblep3.end());
		tempval.it = visiblep3;
                tempval.size = visiblep3.size();
                mapCinsert(t, tempval);

		
	}
	else { 
		parlay::sequence<point> visiblep1;
		t = { p[0],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep1.push_back(p[j]);
			}
		}
		sort(visiblep1.begin(), visiblep1.end());
		tempval.it = visiblep1;
                tempval.size = visiblep1.size();
                mapCinsert(t, tempval);

		parlay::sequence<point> visiblep2;

		t = { p[1],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep2.push_back(p[j]);
			}
		}
		sort(visiblep2.begin(), visiblep2.end());
		tempval.it = visiblep2;
                tempval.size = visiblep2.size();
                mapCinsert(t, tempval);

		parlay::sequence<point> visiblep3;

		t = { p[2],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep3.push_back(p[j]);
			}
		}
		sort(visiblep3.begin(), visiblep3.end());
		tempval.it = visiblep3;
                tempval.size = visiblep3.size();
                mapCinsert(t, tempval);
		
	}

	timer timer1, timer2, timer3, timer4, timer5, timer6, timer7;
	facet tempf;
	//初始M
	timer1.start();
	if (start(p[0], p[1], p[2])) {
		facet11 = { p[1],p[0] };facet12 = { p[0],p[2] };
		facet21 = { p[2],p[1] }; facet22 = { p[1],p[0] };
		facet31 = { p[0],p[2] }; facet32 = { p[2],p[1] };
		cilk_for(int i = 0; i < 3; i++) {
			if (i == 0) {
				ProcessRidge(facet11, p[0], facet12);
			}
			if (i == 1) {
				ProcessRidge(facet21, p[1], facet22);
			}
			if (i == 2) {
				ProcessRidge(facet31, p[2], facet32);
			}			
		}
	}
	else {
		facet11 = { p[2],p[0] }; facet12 = { p[0],p[1] };
		facet21 = { p[0],p[1] }; facet22 = { p[1],p[2] };
		facet31 = { p[1],p[2] }; facet32 = { p[2],p[0] };
		cilk_for(int i = 0; i < 3; i++) {
			if (i == 0) {
				ProcessRidge(facet11, p[0], facet12);
			}
			if (i == 1) {
				ProcessRidge(facet21, p[1], facet22);
			}
			if (i == 2) {
				ProcessRidge(facet31, p[2], facet32);
			}
		}
	}
	timer1.stop();

	int Hsize = 0;
	for (int i = 0; i < n; i++) {
		if (H[i].v1.n != 0) {
			//cout << H[i].v1.pivot << "," << H[i].v2.pivot << endl;
			Hsize++;
		}
		
	}
	/*for (int i = 0; i < n; i++) {
		cout << H[i].v1.pivot << "," << H[i].v2.pivot << endl;
		cout << H[i].v1.n << "," << H[i].v2.n << endl;
	}*/
	cout << "total time: " << timer1.get_total() << endl;
	cout << "output size: " << Hsize << endl;
	cout << "1: " << count1 << endl;
	cout << "2: " << count2 << endl;
	cout << "3: " << count3 << endl;
	cout << "4: " << count4 << endl;
	cout << "5: " << count5 << endl;
	//cout << "time: " << testtime1.get_total() << endl;
	cout << endl;
	//cout << "M size:           			" << M.size() << endl;

	cout << "init and find t1,t2:			" << testtime1.get_total() << endl;
	
	
	cout << "check visible:				" << testtime2.get_total() << endl;
	
	cout << "mapC insert and H insert and erase:	" << testtime3.get_total() << endl;
	
	cout << "M map check:				" << testtime4.get_total() << endl;
	cout << "case 4 total time:			" << testtime5.get_total() << endl;

	/*for (unordered_set<Facet>::iterator it = H.begin(); it != H.end(); ++it) {
		cout << it->v1.pivot << "," << it->v2.pivot << endl;
	}
	*/
	return 0;
}


