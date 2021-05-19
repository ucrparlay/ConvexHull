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
#include "sequence.h"

using namespace std;
using namespace pbbs;

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

int hash3(uint32_t a, uint32_t b, int size) {
	a = (hash2(a) ^ hash2(b) << 1)% size;
	if (a < 0) a = -a;
	return (a);
}

struct point
{
	double x, y;
	int pivot, n;
	bool operator < (const point &o) const
	{
		if (n == 0) return false;
		else if (o.n == 0) return true;
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
	point min;
	int size,flag;
};

facet* H;
point* p;
point* space;
bool* small;
bool* mid;
bool* large;
pair <int, int> *R;
pair <facet, mapCvalue> *hashC;
int* hashlock;
int n10,n100,n1000,n10000,n100000,n1000000;

int hashfacet(facet o) {
	//cout << o.v1.x << "," << o.v1.y << "," << o.v1.n << endl;
	int a = (hash1(o.v1.pivot) ^ hash1(o.v2.pivot) << 1) % (o.v1.n * 10);
	if (a < 0) a = -a;
	return (a);
}


bool mapCinsert(facet t, mapCvalue val) {
	int i = hashfacet(t);//cout << "insert:" << t.v1.pivot << "," << t.v2.pivot << ";" << t.v1.n << endl;
	while (!pbbs::atomic_compare_and_swap(&hashlock[i], 0, 1)) {
		i = (i + 1) % (t.v1.n * 10);
	}
	hashC[i] = pair<facet, mapCvalue>(t, val);
	return true;
}

mapCvalue mapCfind(facet t) {
	int i = hashfacet(t);//cout << "need find:" << t.v1.pivot << "," << t.v2.pivot << endl;
	while (!(hashC[i].first == t)) {//cout << hashC[i].first.v1.pivot << "," << hashC[i].first.v2.pivot << endl;
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

void PercentageRandmOnCircle(point *p, int n, int percent) {
	for (int i = 0; i < n*percent/10; i++) {
                p[i].x = (double)n + n * cos((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
                p[i].y = (double)n + n * sin((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
                p[i].pivot = i;
                p[i].n = n;
        }
        for (int i = n*percent/10; i < n; i++) {
                int temp = hash2(i)%(n-2) + 1 ;
		p[i].x = (double)temp + temp * cos((((hash1(i)) % (temp * 100))*360.0) / (100 * temp)*PI / 180);
                p[i].y = (double)temp + temp * sin((((hash1(i)) % (temp * 100))*360.0) / (100 * temp)*PI / 180);
                p[i].pivot = i;
                p[i].n = n;
        }
	double tempx, tempy;
	for (int i = 0; i < n; i++) {
		tempx = p[hash2(i)%n].x;
		p[hash2(i)%n].x = p[i].x;
		p[i].x = tempx;
		tempy = p[hash2(i)%n].y;
                p[hash2(i)%n].y = p[i].y;
                p[i].y = tempy;
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

int computevisblepoint1(point* A, facet t,int n) {
	if (n <= 256) {
		int ret = 0;
		for (int i = 0; i < n; i++) {
			if (visible(A[i], t)) ret++;
		}
		return ret;
	}
	int L, R;
	L = cilk_spawn computevisblepoint1(A, t, n / 2);
	R = computevisblepoint1(A + n / 2, t, n - n / 2);
	cilk_sync;
	return L + R;
}



void pointmerge(mapCvalue v1, mapCvalue v2, facet t) {
	int pointer0 = 0;
	if (v1.size == 0) {
		pointer0 = computevisblepoint1(v2.it,t,v2.size);
	}
	else if (v2.size == 0) {
		pointer0 = computevisblepoint1(v1.it, t, v1.size);
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
	mapCvalue tempval;
	tempval.it = pmerge;
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
	
	tempval.size = pointer0;
	mapCinsert(t, tempval);
	return;
}

void scan_inplace__ (int* in, int n) {
    if (n <= 10000) {
        for (int i = 1;i < n;i++) {
            in[i] += in[i-1];
        }
        return;
    }
    int root_n = (int)sqrt(n);
    int* offset = new int[root_n-1];    
    cilk_for (int i = 0;i < root_n-1;i++) {
        offset[i] = 0;
        for (int j = i*root_n;j < (i+1)*root_n;j++) {
            offset[i] += in[j];
        }
    }
    for (int i = 1;i < root_n-1;i++) offset[i] += offset[i-1];

    cilk_for (int i = 0;i < root_n;i++) {
        if (i == root_n-1) {
            for (int j = i*root_n+1;j < n;j++) {
                in[j] += in[j-1];
            }
        } else {
            for (int j = i*root_n+1;j < (i+1)*root_n;j++) {
                in[j] += in[j-1];
            }
        }
    }
        
    cilk_for (int i = 1;i < root_n;i++) {
        if (i == root_n-1)  {
            for (int j = i * root_n;j < n;j++){
                in[j] += offset[i-1];
            }
        } else {
            for (int j = i * root_n;j < (i+1)*root_n;j++) {
                in[j] +=  offset[i-1];
            }
        }
    }
    delete[] offset;
}

void prefixsum(int* In, int* Out, int* B, int* C, int n) {
	if (n == 0) return;
	if (n == 1) {
		Out[0] = In[0];
		return;
	}
	cilk_for(int i = 0; i < n / 2; i++)
		B[i] = In[2 * i] + In[2 * i + 1];

	prefixsum(B, C, B + n / 2, C + n / 2, n / 2);
	Out[0] = In[0];

	cilk_for(int i = 1; i < n; i++) {
		if (i % 2) Out[i] = C[i / 2];
		else Out[i] = C[i / 2 - 1] + In[i];
	}
}

int mysum(int size, int* arr) {
	int output = 0;
	for (int i = 0; i < size; i++) {
		output += arr[i];
	}
	return output;
}

void pcopy(point *s, int size, point *news){
	cilk_for (int i = 0; i < size; i++) {
		news[i] = s[i];//cout << i << ":" << s[i].pivot << " ";
	}//cout << endl;
}

void findsmallest(point *arr, int size) {
	for (int i = 1; i < size; i++) {//cout << "n:" << arr[0].n << " pivot:" << arr[0].pivot << endl;
		if(arr[0].n == 0) arr[0] = arr[i];
		else if (arr[i] < arr[0] && arr[i].n != 0) {
			arr[0] = arr[i];
		}
	}
}

timer tim1,tim2,tim3,tim4;

void parallelmerge(mapCvalue v1, mapCvalue v2, facet t, facet t2) {
	int n1,n2,n3,offset = 1000;
	//tim1.start();
	/*if(v1.flag && v2.flag){
		pointmerge(v1, v2, t)
	}*/
	if(v1.size == 0) n1 = 0;
	else n1 = (v1.size-1) / offset + 1;
	//cout << "ms:" << t.v1.pivot << ";" << t.v2.pivot << endl;
	if(v2.size == 0) n2 = 0;
	else n2 = (v2.size-1) / offset + 1;
	/*if((v1.size+v2.size)<100) n100++;
	else if((v1.size+v2.size)<1000) n1000++;
        else if((v1.size+v2.size)<10000) n10000++;
        else if((v1.size+v2.size)<100000) n100000++;
	else n1000000++;*/

	//cout << "n1,n2:" << n1 << "," << n2 << endl;
	//sequence<int> count1(n1,0);
	//sequence<int> count2(n2,0);
	int *count1 = new int[n1];
	int *count2 = new int[n2];
/*	int *count1o = new int[n1];
	int *count2o = new int[n2];
	int *count1b = new int[n1];
	int *count2b = new int[n2];
	int *count1c = new int[n1];
	int *count2c = new int[n2];*/
	//point zeropoint = {0,0,0,0};
	//sequence<point> smallest(n1+n2);
	point *smallest = new point[n1+n2];
	//tim1.stop();
	point zeropoint = {0,0,0,0};
	//cout << "init smallest:";
	//tim2.start();
        cilk_for(int i = 0;i<n1+n2;i++) smallest[i] = zeropoint;
	//tim2.stop();//cout << smallest[i].pivot << ",";cout << endl;
	//tim1.start();
	//point *vpoint1 = new point[v1.size];
	//point *vpoint2 = new point[v2.size];
	//tim1.stop();
	//cout << "v1.size:" << v1.size << ";v2.size:" << v2.size << endl;
	cilk_for (int i = 0; i < n1; i++) {
		int number = 0;
		for (int j = i*offset; j < i*offset + offset && j < v1.size; j++) {
			if (visible(v1.it[j], t) && !visible(v1.it[j], t2)) {
				point temp = v1.it[j];
				if (temp < smallest[i] || number == 0) {
					smallest[i] = temp;
				}
				//vpoint1[i*offset + number] = temp;
				number++;
			}//cout << number << endl;
		}
		count1[i] = number;
	}//cout << "check1" << "," << v1.size << "," << n1 << "," << sizeof(count1) << "," << count1[0] << endl;
	cilk_for (int i = 0; i < n2; i ++) {
		int number = 0;
		for (int j = i*offset; j < i*offset + offset && j < v2.size; j++) {
			if (visible(v2.it[j], t)) {
				point temp = v2.it[j];
				if (temp < smallest[n1+(i)] || number == 0) {
					smallest[n1 + (i)] = temp;
				}
				//vpoint2[i*offset + number] = temp;
				number++;
			}
		}
		count2[i] = number;//cout << "," << count1[0] << endl;
	}//cout << "check2" << "," << count2[0] << endl;
	//cout << sizeof(count1) << "," << count1[0] << "," << count1[1] << endl;
	/*cout << "v1:";
        for(int i = 0;i<v1.size;i++) cout << v1.it[i].pivot << ",";cout << endl;
	cout << "v1:";
        for(int i = 0;i<v2.size;i++) cout << v2.it[i].pivot << ",";cout << endl;
	cout << "all smallest:";
	for(int i = 0;i<n1+n2;i++) cout << smallest[i].pivot << ",";cout << endl;
	*/
	//for(int i = 0;i<n1;i++) cout << " " << count1[i];cout << endl;
	//tim3.start();
	scan_inplace__(count1,n1);//for(int i = 0;i<n1;i++) cout << " " << count1[i];cout << endl;
	scan_inplace__(count2,n2);
	//tim3.stop();
	findsmallest(smallest,n1 + n2);
	int n2size = count2[n2-1];
	int n1size = count1[n1-1];//cout << "check3" << endl;
	if(n1size + n2size == 0){
		mapCvalue tempval;
        	tempval.size = 0;
        	mapCinsert(t, tempval);
        return;
	}
	int psize = n1size + n2size;//cout << psize << endl;
	point *pmerge = new point[psize];//cout << "check5" << endl;
	//point *tempmerge = pmerge;
	//cout << sizeof(count1) << "," << count1[0] << "," << count1[1] << endl;
	cilk_for (int i = 0; i < n1; i++) {
         	int number = 0;     
                for (int j = i*offset; j < i*offset + offset && j < v1.size; j++) {
                        if (visible(v1.it[j], t) && !visible(v1.it[j], t2)) {
                        	if(i == 0)pmerge[number] = v1.it[j];
				else pmerge[count1[i-1] + number] = v1.it[j];
			 	number++;
                        }
                }
               
        }
        cilk_for (int i = 0; i < n2; i ++) {
                int number = 0;
                for (int j = i*offset; j < i*offset + offset && j < v2.size; j++) {
                        if (visible(v2.it[j], t)) {
            			if(i == 0)pmerge[count1[n1-1] + number] = v2.it[j];
                                else pmerge[count1[n1-1] + count2[i-1] + number] = v2.it[j];
                                number++;                    
                        }
                }
             
        }
	/*
	cilk_for (int i = 0; i < n1; i++) {
		//cout << "check6" << endl;
		if (i == 0) {
                        pcopy(vpoint1, count1[0], pmerge);
                }
                else {
                        pcopy(vpoint1 + i * offset, count1[i] - count1[i-1], pmerge + count1[i-1]);
                }

	}//cout << "check4" << endl;
	cilk_for (int i = 0; i < n2; i++) {
		if (i == 0) {
			pcopy(vpoint2, count2[0], pmerge + count1[n1-1]);
		}
		else {
			pcopy(vpoint2 + i * offset, count2[i] - count2[i-1], pmerge + count1[n1-1] + count2[i - 1]);
		}
	}*/
	//cout << "end" << endl;
	//cout << "psize:" << psize << endl;
	//cout << "1:";
	//for(int i = 0;i<n1;i++) cout << count1[i] << ",";cout << endl;cout << "2:";
	//for(int i = 0;i<n2;i++) cout << count2[i] << ",";cout << endl;
	//for(int i = 0;i<psize;i++)cout << pmerge[i].pivot << ",";
	//cout << endl;
	mapCvalue tempval;//cout << "end" << endl;
	tempval.min = smallest[0];//cout << "smallest:" << smallest[0].pivot << endl;
	tempval.it = pmerge;//cout << "end" << endl;
	tempval.size = psize;//cout << "end" << endl;
	tempval.flag = 0;
	mapCinsert(t, tempval);//cout << "end" << endl;
	return;
}

void ProcessRidge(facet t1, point r, facet t2) {
	
	//cout << "r:" << t1.v1.pivot << "," << r.pivot << "," << t2.v2.pivot << endl;
	mapCvalue t1val = mapCfind(t1);
        mapCvalue t2val = mapCfind(t2);
	//cout << t1val.size << "," << t2val.size << endl;
	//cout << "t:";
	//for(int i = 0;i < t2val.size;i++){
	//cout << t2val.it[i].pivot << " ";}
	//cout << endl;
	if ((t1val.size == 0) && (t2val.size == 0)) {
		H[t1.v1.pivot] = t1;
	}
	else if (!(t2val.size == 0) && !(t1val.size == 0) && t1val.min == t2val.min) {
	}
	else if ((t1val.size == 0) || (t2val.size != 0 && (t2val.min < t1val.min))) {
		facet t;
		t.v1 = r;
		t.v2 = t2val.min;
		//cout << "start" << endl;	
		parallelmerge(t2val, t1val, t, t1);//cout << "merge finish" << endl;
		cilk_for(int i = 0; i < 2; i++) {
			if (i) {
				ProcessRidge(t1, r, t);
			}
			else if (!InsertandSetL(t)) {//cout << "checkL" << endl;
				facet oldt = getvalueL(t);//cout << "checkget" << endl;
				ProcessRidge(t, t.v2, oldt);
			}
		}
	}
	else {
		facet t;
		t.v1 = t1val.min;
		t.v2 = r;
		//cout << "start" << endl;
		parallelmerge(t1val, t2val, t, t2);//cout << "merge finish" << endl;
		cilk_for(int i = 0; i < 2; i++) {
			if (i) {
				ProcessRidge(t, r, t2);
			}
			else if (!InsertandSetR(t)) {//cout << "checkR" << endl;
				facet oldt = getvalueR(t);//cout << "checkget" << endl;
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
	else if (in == 10) {
		RandmOnCircle(p, n);
	}
	else if (in == 1) PercentageRandmOnCircle(p,n,1);
        else if (in == 2) PercentageRandmOnCircle(p,n,2);
        else if (in == 3) PercentageRandmOnCircle(p,n,3);
        else if (in == 4) PercentageRandmOnCircle(p,n,4);
        else if (in == 5) PercentageRandmOnCircle(p,n,5);
        else if (in == 6) PercentageRandmOnCircle(p,n,6);
        else if (in == 7) PercentageRandmOnCircle(p,n,7);
        else if (in == 8) PercentageRandmOnCircle(p,n,8);
        else if (in == 9) PercentageRandmOnCircle(p,n,9);

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
		if (!initpointer) {
			tempval.size = initpointer;
			mapCinsert(t, tempval);
		}
		else{
			point *visiblep1 = new point[initpointer];
			initpointer = 0;
			for (int j = 3; j < n; j++) {
				if (visible(p[j], t)) {
					visiblep1[initpointer] = p[j];
					initpointer++;
				}
			}
			tempval.min = visiblep1[0];
			tempval.it = visiblep1;
			tempval.size = initpointer;
			tempval.flag = 0;
			mapCinsert(t, tempval); 
		}

		initpointer = 0;
		t = { p[2],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		if (!initpointer) {
			tempval.size = initpointer;
			mapCinsert(t, tempval);
		}
		else {
			point *visiblep2 = new point[initpointer];
			initpointer = 0;
			for (int j = 3; j < n; j++) {
				if (visible(p[j], t)) {
					visiblep2[initpointer] = p[j];
					initpointer++;
				}
			}
			tempval.min = visiblep2[0];
			tempval.it = visiblep2;
			tempval.size = initpointer;
                        tempval.flag = 0;
			mapCinsert(t, tempval);
		}

		initpointer = 0;
		t = { p[1],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		if (!initpointer) {
			tempval.size = initpointer;
			mapCinsert(t, tempval);
		}
		else {
			point *visiblep3 = new point[initpointer];
			initpointer = 0;
			for (int j = 3; j < n; j++) {
				if (visible(p[j], t)) {
					visiblep3[initpointer] = p[j];
					initpointer++;
				}
			}
			tempval.min = visiblep3[0];
			tempval.it = visiblep3;
			tempval.size = initpointer;
                        tempval.flag = 0;
			mapCinsert(t, tempval);
		}
	}
	else {
		int initpointer = 0;
		t = { p[0],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		if (!initpointer) {
			tempval.size = initpointer;
			mapCinsert(t, tempval);
		}
		else {
			point *visiblep1 = new point[initpointer];
			initpointer = 0;
			for (int j = 3; j < n; j++) {
				if (visible(p[j], t)) {
					visiblep1[initpointer] = p[j];
					initpointer++;
				}
			}
			tempval.min = visiblep1[0];
			tempval.it = visiblep1;
			tempval.size = initpointer;
                        tempval.flag = 0;
			mapCinsert(t, tempval);
		}

		initpointer = 0;
		t = { p[1],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		if (!initpointer) {
			tempval.size = initpointer;
			mapCinsert(t, tempval);
		}
		else {
			point *visiblep2 = new point[initpointer];
			initpointer = 0;
			for (int j = 3; j < n; j++) {
				if (visible(p[j], t)) {
					visiblep2[initpointer] = p[j];
					initpointer++;
				}
			}
			tempval.min = visiblep2[0];
			tempval.it = visiblep2;
			tempval.size = initpointer;
                        tempval.flag = 0;
			mapCinsert(t, tempval);
		}

		initpointer = 0;
		t = { p[2],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				initpointer++;
			}
		}
		if (!initpointer) {
			tempval.size = initpointer;
			mapCinsert(t, tempval);
		}
		else {
			point *visiblep3 = new point[initpointer];
			initpointer = 0;
			for (int j = 3; j < n; j++) {
				if (visible(p[j], t)) {
					visiblep3[initpointer] = p[j];
					initpointer++;
				}
			}
			tempval.min = visiblep3[0];
			tempval.it = visiblep3;
			tempval.size = initpointer;
                        tempval.flag = 0;
			mapCinsert(t, tempval);
		}
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
	space = new point[96*n];
	small = new bool[32 * n];
	mid = new bool[32 * n / 100];
	large = new bool[32 * n / 10000];
	R = new pair<int, int>[10 * n];
	H = new facet[n];
	hashC = new pair<facet, mapCvalue>[10 * n];
	hashlock = new int[10 * n];
	timer timer1;
	facet facet11, facet12, facet21, facet22, facet31, facet32;
	initpoints(type_of_input, n);cout << "check" << endl;
	timer1.start();
	
	init(n);//cout << "check" << endl;
	//for(int i = 0; i < n; i++){cout << p[i].x << "," << p[i].y << endl;}
	if (start(p[0], p[1], p[2])) {
		//cout << "check" << endl;
		facet11 = { p[1],p[0] }; facet12 = { p[0],p[2] };
		facet21 = { p[2],p[1] }; facet22 = { p[1],p[0] };
		facet31 = { p[0],p[2] }; facet32 = { p[2],p[1] };//cout << "check" << endl;
		cilk_for(int i = 0; i < 3; i++) {
			if (i == 0) {//cout << "check" << endl;
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
	else {//cout << "check" << endl;
		facet11 = { p[2],p[0] }; facet12 = { p[0],p[1] };
		facet21 = { p[0],p[1] }; facet22 = { p[1],p[2] };
		facet31 = { p[1],p[2] }; facet32 = { p[2],p[0] };//cout << "check" << endl;
		cilk_for(int i = 0; i < 3; i++) {
			if (i == 0) {//cout << "check" << endl;
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

	//cout << "t1:" << tim1.get_total() << "\nt2:" << tim2.get_total() << "\nt3:" << tim3.get_total() << "\nt4:" << tim4.get_total() << endl;
	//cout << "100:" << n100 << "\n1000:" << n1000 << "\n10000:" << n10000 << "\n100000:" << n100000 << "\n1000000:" << n1000000 << endl;
	cout << "total time: " << timer1.get_total() << endl;
	cout << "output size: " << Hsize << endl;
	return 0;
}
