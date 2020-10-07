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

struct parallelmap
{
	
};

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
	parlay::sequence<point>::iterator it;
	int size;
};

int hashfacet(facet o) {
	return (hash1(o.v1.pivot) ^ (hash1(o.v2.pivot) << 1) % (o.v1.n*10));
}

/*namespace std {
	template<>
	struct hash<Facet> {
		std::size_t operator()(const Facet& o) const{
			using std::hash;
			return (hash<int>{}(o.v1.pivot) ^ (hash<int>{}(o.v2.pivot) << 1));
		}
	};
	template<>
	struct hash<point> {
		std::size_t operator()(const point& o) const {
			using std::hash;
			return (hash<int>{}(o.pivot));
		}
	};
}*/
using namespace std;

//unordered_map<facet, parlay::sequence<point>> mapC;
//unordered_map<point, facet> M;
facet* H;
point* p;
pair <int, int>* R;
pair<facet,mapCvalue>* hashC;
parlay::hashtable<parlay::hash_numeric<int>> table(10000000, parlay::hash_numeric<int>{});
//point zerop;
//facet zerof;


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
	//cout << "check ProcessRidge 1" << endl;
	testtime1.start();
	//unordered_map<facet, parlay::sequence<point>>::iterator t1iter;
	//unordered_map<facet, parlay::sequence<point>>::iterator t2iter;
	parlay::sequence<point>::iterator t1iter;
	parlay::sequence<point>::iterator t2iter;
	parlay::sequence<point> vp;
	int t1size, t2size;
	//cout << t1.v1.pivot << "," << t1.v2.pivot << " " << t2.v1.pivot << "," << t2.v2.pivot << endl;
	mapCvalue t1val = mapCfind(t1);
	t1iter = t1val.it;
	mapCvalue t2val = mapCfind(t2);
	t2iter = t2val.it;
	//cout << "check ProcessRidge 2" << endl;
	testtime1.stop();
	if ((t1val.size == 0) && (t2val.size == 0)) {
		//cout << "1" << endl;
		//cout << t1.v1.pivot << "," << t1.v2.pivot << " " << t2.v1.pivot << "," << t2.v2.pivot << endl;
		H[t1.v1.pivot] = t1;
		count1++;
		
	}
	else if (!(t2val.size == 0) &&  !(t1val.size == 0) && *(t1val.it) == *(t2val.it)) {
		//cout << "2" << endl;
		count2++;
		testtime1.start();
		//H[t1.v1.pivot] = zerof;
		//H[t2.v1.pivot] = zerof;
		testtime1.stop();
	}
	else if ((t1val.size == 0) || (t2val.size != 0 && (*(t2val.it) < *(t1val.it)))) {//左
		//cout << "3" << endl;
		testtime2.start();
		count3++;
		parlay::sequence<point> c;
		facet t;
		t.v1 = r;
		t.v2 = *(t2val.it);
		//cout << "3-2: " << t2iter->second.size() << endl;
		//cout << "3-1: " << t1iter->second.size() << endl;
		if (t1val.size == 0) {
			//sort(t2iter->second.begin(), t2iter->second.end());
			for (int i = 0; i < t2val.size; i++) {
				//cout << i << endl;
				count5++;
				if (visible(*t2iter, t)) {
					vp.push_back((*t2iter);
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
					if (visible(*(t1iter + pointer1), t)) {
						vp.push_back(*(t1iter + pointer1));
					}
				}
			}
			else if (pointer2 != t2val.size) {
				for (int i = pointer2; i < t2val.size; i++) {
					if (visible(*(t2iter + pointer2), t)) {
						vp.push_back(*(t2iter + pointer2));
					}
				}
			}
			//sort(t1iter->second.begin(), t1iter->second.end());
			//sort(t2iter->second.begin(), t2iter->second.end());
			/*set_union(t2iter->second.begin(), t2iter->second.end(), t1iter->second.begin(), t1iter->second.end(), inserter(c, c.begin()));
			for (int i = 0; i < c.size(); i++) {
				//cout << i << endl;
				count5++;
				if (visible(c[i], t)) {
					vp.push_back(c[i]);
				}
			}*/
		}
		/*for (int i = 0; i < t2iter->second.size(); i++) {
			//cout << i << endl;
			count5++;
			if (visible(t2iter->second[i], t)) {
				if (vp.empty()) {
					vp.push_back(t2iter->second[i]);
				}
				else {
					if (t2iter->second[i] < vp.back()) {
						vp.push_back(t2iter->second[i]);
					}
					else {
						temp = vp.back();
						vp[vp.size() - 1] = t2iter->second[i];
						vp.push_back(temp);
					}
				}
			}
		}
		
		for (int i = 0; i < t1iter->second.size(); i++) {
			//cout << i << endl;
			count5++;
			if (visible(t1iter->second[i], t)) {
				if (vp.empty()) {
					vp.push_back(t1iter->second[i]);
				}
				else {
					if (t1iter->second[i] < vp.back()) {
						vp.push_back(t1iter->second[i]);
					}
					else {
						temp = vp.back();
						vp[vp.size() - 1] = t1iter->second[i];
						vp.push_back(temp);
					}
				}
			}
		}*/
		
		testtime2.stop();
		testtime3.start();
		mapCinsert(t, vp);

		
		//H[t2.v1.pivot] = zerof;
		//H[t.v1.pivot] = t;

		//cout << t.v1.pivot << "," << t.v2.pivot << endl;
		testtime3.stop();
		for (int i = 0; i < 2; i++) {
			//cout << "check3:" << i << endl;
			if (i) {
				ProcessRidge(t1, r, t);
			}
			else if(!InsertandSetL(t)){
				//cout << "check3 1" << endl;
				facet oldt = getvalueL(t);

				//cout << "check3 1" << endl;
				ProcessRidge(t, t.v2, oldt);
			}
		}
	}
	else {//右
		//cout << "4" << endl;
		testtime5.start();
		count4++;
		parlay::sequence<point> c;
		facet t;
		t.v1 = *(t1val.it);
		t.v2 = r;
		//cout << "4-2: " << t2iter->second.size() << endl;
		//cout << "4-1: " << t1iter->second.size() << endl;
		if (t2iter->second.empty()) {
			//sort(t1iter->second.begin(), t1iter->second.end());
			for (int i = 0; i < t1val.size; i++) {
				//cout << i << endl;
				count5++;
				if (visible(*(t1iter + i), t)) {
					vp.push_back(*(t1iter + i));
				}
			}
		}
		else {
			int pointer1 = 0;
			int pointer2 = 0;
			while (pointer1 < t1val.size && t2val.size)
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
			//sort(t1iter->second.begin(), t1iter->second.end());
			//sort(t2iter->second.begin(), t2iter->second.end());
			/*set_union(t2iter->second.begin(), t2iter->second.end(), t1iter->second.begin(), t1iter->second.end(), inserter(c, c.begin()));
			for (int i = 0; i < c.size(); i++) {
				//cout << i << endl;
				count5++;
				if (visible(c[i], t)) {
					vp.push_back(c[i]);
				}
			}*/
		}
		
		/*for (int i = 0; i < t2iter->second.size(); i++) {
			//cout << i << endl;
			count5++;
			if (visible(t2iter->second[i], t)) {
				if (vp.empty()) {
					vp.push_back(t2iter->second[i]);
				}
				else {
					if (t2iter->second[i] < vp.back()) {
						vp.push_back(t2iter->second[i]);
					}
					else {
						temp = vp.back();
						vp[vp.size() - 1] = t2iter->second[i];
						vp.push_back(temp);
					}
				}
			}
		}
		
		for (int i = 0; i < t1iter->second.size(); i++) {
			//cout << i << endl;
			count5++;
			if (visible(t1iter->second[i], t)) {
				if (vp.empty()) {
					vp.push_back(t1iter->second[i]);
				}
				else {
					if (t1iter->second[i] < vp.back()) {
						vp.push_back(t1iter->second[i]);
					}
					else {
						temp = vp.back();
						vp[vp.size() - 1] = t1iter->second[i];
						vp.push_back(temp);
					}
				}
			}
		}*/
		//cout << "4check4" << endl;
		mapCinsert(t, vp);

		//H[t1.v1.pivot] = zerof;
		//H[t.v1.pivot] = t;
		
		point tempr;
		facet tempt;
		testtime5.stop();
		for (int i = 0; i < 2; i++) {
			//cout << "check4:" << i << endl;
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
	parlay::sequence<point>::iterator it;
	cout << &hashtest << endl;
	cout << sizeof(it) << endl;
	/*
	zerop.x = 0;
	zerop.y = 0;
	zerop.pivot = 0;
	zerop.n = 0;
	zerof.v1 = zerop;
	zerof.v2 = zerop;
	*/
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

	/*for (int i = 0; i < n*100; i++) {
		cout << R[i].first << "," << R[i].second << endl;
	}
	for (int i = 0; i < n; i++) {
		cout << p[i].x << "," << p[i].y << endl;
	}*/


	while (collinear(p[0], p[1], p[2])) {//排除初始三点共线
		int x = rand() % n;
		point temp = p[2];
		p[2] = p[x];
		p[x] = temp;
	}
	parlay::sequence<point> visiblep;
	point temp;
	facet t, facet11, facet12, facet21, facet22, facet31, facet32;
	//确定每条线中的点顺时针排列
	//
	//cout << "初始C,H" << endl;
	if (start(p[0], p[1], p[2])) {
		t = { p[0],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.push_back(p[j]);
			}
		}
		sort(visiblep.begin(), visiblep.end());
		mapCinsert(t, visiblep);
		//H[t.v1.pivot] = t;
		visiblep.clear();
		t = { p[2],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.push_back(p[j]);
			}
		}
		sort(visiblep.begin(), visiblep.end());
		mapCinsert(t, visiblep);
		//H[t.v1.pivot] = t;
		visiblep.clear();
		t = { p[1],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.push_back(p[j]);
			}
		}
		sort(visiblep.begin(), visiblep.end());
		mapCinsert(t, visiblep);
		//H[t.v1.pivot] = t;
		visiblep.clear();
	}
	else {
		t = { p[0],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.push_back(p[j]);
			}
		}
		sort(visiblep.begin(), visiblep.end());
		mapCinsert(t, visiblep);
		//H[t.v1.pivot] = t;
		visiblep.clear();

		t = { p[1],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.push_back(p[j]);
			}
		}
		sort(visiblep.begin(), visiblep.end());
		mapCinsert(t, visiblep);
		//H[t.v1.pivot] = t;
		visiblep.clear();

		t = { p[2],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.push_back(p[j]);
			}
		}
		sort(visiblep.begin(), visiblep.end());
		mapCinsert(t, visiblep);
		//H[t.v1.pivot] = t;
		visiblep.clear();
	}

	timer timer1, timer2, timer3, timer4, timer5, timer6, timer7;
	facet tempf;
	//初始M
	/*for (int i = 0; i < 3; i++) {

		if (start(p[0], p[1], p[2])) {
			tempf = {p[i] ,p[(i+2)%3] };
			M.insert(pair<point, facet>(p[i], tempf));

			tempf = { p[(i+1)%3],p[i] };
			M.insert(pair<point, facet>(p[i], tempf));
		}
		else {
			tempf = { p[(i + 2) % 3], p[i] };
			M.insert(pair<point, facet>(p[i], tempf));

			tempf = { p[i], p[(i + 1) % 3]};
			M.insert(pair<point, facet>(p[i], tempf));
		}
	}*/
	timer1.start();
	if (start(p[0], p[1], p[2])) {
		facet11 = { p[1],p[0] };facet12 = { p[0],p[2] };
		facet21 = { p[2],p[1] }; facet22 = { p[1],p[0] };
		facet31 = { p[0],p[2] }; facet32 = { p[2],p[1] };
		for(int i = 0; i < 3; i++) {
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
		for(int i = 0; i < 3; i++) {
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
	cout << "mapc size: " << mapC.size() << endl;
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


