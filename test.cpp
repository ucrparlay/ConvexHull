#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <atomic>
#include <mutex>
#include <string.h>
#include <cmath>
#include <algorithm>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include "get_time.h"
#include "utilities.h"


using namespace std;

const  double PI = 3.1415926;
const  double EPSINON = 0.0000001;

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

inline uint32_t hash2(uint32_t a) {
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
	int pivot, n;
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
	point *it;
	int size;
};

facet* H;
point* p;
pair <int, int> *R;
pair <facet, mapCvalue> *hashC;
int* hashlock;
//point* pmerge;

int hashfacet(facet o) {
	int a = (hash1(o.v1.pivot) ^ hash1(o.v2.pivot) << 1) % (o.v1.n * 10);
	if (a < 0) a = -a;
	return (a);
}


bool mapCinsert(facet t, mapCvalue val) {
	int i = hashfacet(t);
	while (!pbbs::atomic_compare_and_swap(&hashlock[i], 0, 1)) {
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
	while (!pbbs::atomic_compare_and_swap(&R[i], pair <int, int>(0, 0), pair<int, int>(t.v1.pivot, t.v2.pivot))) {
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

facet getvalueL(facet t) {
	int i = hash1(t.v2.pivot) % (t.v2.n * 10);
	while (!(R[i].first == t.v2.pivot)) {
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
		p[i].x = (double)n + n * cos((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
		p[i].y = (double)n + n * sin((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
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
	if ((v1.x - v3.x) * (v2.y - v3.y) - (v1.y - v3.y) * (v2.x - v3.x) < EPSINON && (v1.x - v3.x) * (v2.y - v3.y) - (v1.y - v3.y) * (v2.x - v3.x) > -EPSINON) return true;
	else return false;
}

void pointmerge(mapCvalue v1, mapCvalue v2, facet t) {
	int pointer0 = 0;
	if (v1.size == 0) {
		for (int i = 0; i < v2.size; i++) {
			if (visible(v2.it[i], t)) {
				pointer0++;
			}
		}
	}
	else if (v2.size == 0) {
		for (int i = 0; i < v1.size; i++) {
			if (visible(v1.it[i], t)) {
				pointer0++;
			}
		}
	}
	else {
		int pointer1 = 0;
		int pointer2 = 0;
		while (pointer1 < v1.size && pointer2 < v2.size)
		{
			if (v1.it[pointer1] < v2.it[pointer2]) {
				if (visible(v1.it[pointer1], t)) {
					pointer0++;
				}
				pointer1++;
			}
			else if (v2.it[pointer2] < v1.it[pointer1]) {
				if (visible(v2.it[pointer2], t)) {
					pointer0++;
				}
				pointer2++;
			}
			else {
				if (visible(v1.it[pointer1], t)) {
					pointer0++;
				}
				pointer1++;
				pointer2++;
			}
		}
		if (pointer1 != v1.size) {
			for (int i = pointer1; i < v1.size; i++) {
				if (visible(v1.it[i], t)) {
					pointer0++;
				}
			}
		}
		else if (pointer2 != v2.size) {
			for (int i = pointer2; i < v2.size; i++) {
				if (visible(v2.it[i], t)) {
					pointer0++;
				}
			}
		}
	}
	point *pmerge = new point[pointer0];
	pointer0 = 0;
	if (v1.size == 0) {
		for (int i = 0; i < v2.size; i++) {
			if (visible(v2.it[i], t)) {
				pmerge[pointer0] = v2.it[i];
				pointer0++;
			}
		}
	}
	else if (v2.size == 0) {
		for (int i = 0; i < v1.size; i++) {
			if (visible(v1.it[i], t)) {
				pmerge[pointer0] = v1.it[i];
				pointer0++;
			}
		}
	}
	else {
		int pointer1 = 0;
		int pointer2 = 0;
		while (pointer1 < v1.size && pointer2 < v2.size)
		{
			if (v1.it[pointer1] < v2.it[pointer2]) {
				if (visible(v1.it[pointer1], t)) {
					pmerge[pointer0] = v1.it[pointer1];
					pointer0++;
				}
				pointer1++;
			}
			else if (v2.it[pointer2] < v1.it[pointer1]) {
				if (visible(v2.it[pointer2], t)) {
					pmerge[pointer0] = v2.it[pointer2];
					pointer0++;
				}
				pointer2++;
			}
			else {
				if (visible(v1.it[pointer1], t)) {
					pmerge[pointer0] = v1.it[pointer1];
					pointer0++;
				}
				pointer1++;
				pointer2++;
			}
		}
		if (pointer1 != v1.size) {
			for (int i = pointer1; i < v1.size; i++) {
				if (visible(v1.it[i], t)) {
					pmerge[pointer0] = v1.it[i];
					pointer0++;
				}
			}
		}
		else if (pointer2 != v2.size) {
			for (int i = pointer2; i < v2.size; i++) {
				if (visible(v2.it[i], t)) {
					pmerge[pointer0] = v2.it[i];
					pointer0++;
				}
			}
		}
	}
	mapCvalue tempval;
	tempval.it = pmerge;
	tempval.size = pointer0;
	mapCinsert(t, tempval);
	return;
}

void ProcessRidge(facet t1, point r, facet t2) {
	mapCvalue t1val = mapCfind(t1);
	mapCvalue t2val = mapCfind(t2);
	if ((t1val.size == 0) && (t2val.size == 0)) {
		H[t1.v1.pivot] = t1;
	}
	else if (!(t2val.size == 0) && !(t1val.size == 0) && t1val.it[0] == t2val.it[0]) {
	}
	else if ((t1val.size == 0) || (t2val.size != 0 && (t2val.it[0] < t1val.it[0]))) {//左
		facet t;
		t.v1 = r;
		t.v2 = t2val.it[0];
		pointmerge(t1val, t2val, t);
		cilk_for(int i = 0; i < 2; i++) {
			if (i) {
				ProcessRidge(t1, r, t);
			}
			else if (!InsertandSetL(t)) {
				facet oldt = getvalueL(t);
				ProcessRidge(t, t.v2, oldt);
			}
		}
	}
	else {//右
		facet t;
		t.v1 = t1val.it[0];
		t.v2 = r;
		pointmerge(t1val, t2val, t);
		cilk_for(int i = 0; i < 2; i++) {
			if (i) {
				ProcessRidge(t, r, t2);
			}
			else if (!InsertandSetR(t)) {
				facet oldt = getvalueR(t);
				ProcessRidge(oldt, t.v1, t);
			}
		}
	}
	return;
}

void initpoints(int in, int n) {
	if (in == 0) {
		for (int i = 0; i < n; i++) {
			p[i].x = (hash1(i)) % (n * 2);
			p[i].y = (hash2(i)) % (n * 2);
			p[i].pivot = i;
			p[i].n = n;
		}
	}
	else if (in == 1) {
		RandmOnCircle(p, n);
	}
	else 
		cout << "Usage: ./qsort [num_elements] [0 or 1]" << endl;

}

void init(int n) {
	mapCvalue tempval;
	facet t;
	if (start(p[0], p[1], p[2])) {
		int initpointer = 0;
		t = { p[0],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		point *visiblep1 = new point[initpointer];
		initpointer = 0;
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep1[initpointer] = p[j];
				initpointer++;
			}
		}
		tempval.it = visiblep1;
		tempval.size = initpointer;
		mapCinsert(t, tempval);

		initpointer = 0;
		t = { p[2],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		point *visiblep2 = new point[initpointer];
		initpointer = 0;
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep2[initpointer] = p[j];
				initpointer++;
			}
		}
		tempval.it = visiblep2;
		tempval.size = initpointer;
		mapCinsert(t, tempval);

		initpointer = 0;
		t = { p[1],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		point *visiblep3 = new point[initpointer];
		initpointer = 0;
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep3[initpointer] = p[j];
				initpointer++;
			}
		}
		tempval.it = visiblep3;
		tempval.size = initpointer;
		mapCinsert(t, tempval);

	}
	else {
		int initpointer = 0;
		t = { p[0],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		point *visiblep1 = new point[initpointer];
		initpointer = 0;
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep1[initpointer] = p[j];
				initpointer++;
			}
		}
		tempval.it = visiblep1;
		tempval.size = initpointer;
		mapCinsert(t, tempval);

		initpointer = 0;
		t = { p[1],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		point *visiblep2 = new point[initpointer];
		initpointer = 0;
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep2[initpointer] = p[j];
				initpointer++;
			}
		}
		tempval.it = visiblep2;
		tempval.size = initpointer;
		mapCinsert(t, tempval);

		initpointer = 0;
		t = { p[2],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		point *visiblep3 = new point[initpointer];
		initpointer = 0;
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep3[initpointer] = p[j];
				initpointer++;
			}
		}
		tempval.it = visiblep3;
		tempval.size = initpointer;
		mapCinsert(t, tempval);
	}
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
	hashlock = new int[10 * n];
	timer timer1;
	facet facet11, facet12, facet21, facet22, facet31, facet32;
	initpoints(type_of_input, n);
	timer1.start();
	
	init(n);
	if (start(p[0], p[1], p[2])) {
		facet11 = { p[1],p[0] }; facet12 = { p[0],p[2] };
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
			cout << H[i].v1.pivot << "," << H[i].v2.pivot << endl;
			Hsize++;
		}

	}

	cout << "total time: " << timer1.get_total() << endl;
	cout << "output size: " << Hsize << endl;
	return 0;
}
