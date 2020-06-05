#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <list>
#include <map>
//#include <hash_map>
#include <set>
#include <atomic>
#include <mutex>
#include <string.h>
#include <cmath>
#include <algorithm>
#include "get_time.h"

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

inline uint32_t hash2(uint32_t a)
{
	a = (a + 0x1a976ce5) + (a << 15);
	a = (a ^ 0x98bf7c21) ^ (a >> 13);
	a = (a + 0xc092738a) + (a << 7);
	a = (a + 0xbbc2976e) ^ (a << 9);
	a = (a + 0xf6efd123) + (a << 11);
	a = (a ^ 0x987abc01) ^ (a >> 2);
	return a;
}

//bool *C = new bool[1000001]();//0000


struct point
{
	double x, y;
	int pivot;
	bool operator < (const point &o) const
	{
		return pivot < o.pivot;
	}
	/*bool operator == (const point &o) const
	{
		return pivot == o.pivot;
	}*/
}s[10000000];//

typedef struct facet
{
	point v1, v2;
	bool operator < (const facet &o) const
	{
		return v1.pivot < o.v1.pivot;
	}
}Facet;

struct pairfacet
{
	facet f1, f2;
};
struct  ridgefacet
{
	facet f;
	point r;

};

vector<Facet> facet2;
map<facet, set<point>> mapC;
map<point, set<facet>> InversemapC;
multimap<point, facet> M;
set<Facet> H;
//hash_map<int, ridgefacet> R;

void RandmOnCircle(point *p, int n) {
	for (int i = 0; i < n; i++) {
		p[i].x = (double)n + n * cos((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
		p[i].y = (double)n + n * sin((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
		p[i].pivot = i;
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

double cross_product(point a, point b, point c) {
	return (b.x - a.x)*(c.y - a.y) - (c.x - a.x)*(b.y - a.y);
}

bool cmp1(point a, point b) {
	if (a.y == b.y)
		return a.x < b.x;
	return a.y < b.y;
}

bool cmp2(point a, point b) {
	if (atan2(a.y - s[0].y, a.x - s[0].x) != atan2(b.y - s[0].y, b.x - s[0].x))
		return (atan2(b.y - s[0].y, b.x - s[0].x) - (atan2(a.y - s[0].y, a.x - s[0].x)) > EPSINON);
	return (b.x - a.x > EPSINON);
}

int main(int argc, char** argv) {

	if (argc != 3) {
		cout << "Usage: ./qsort [num_elements] [type of input]" << endl;
		return 0;
	}
	int n = atoi(argv[1]);
	int type_of_input = atoi(argv[2]);
	point* p = new point[n];

	if (type_of_input == 0) {
		for (int i = 0; i < n; i++) {
			p[i].x = (hash1(i)) % (n * 2);
			p[i].y = (hash2(i)) % (n * 2);
			p[i].pivot = i;
		}
	}
	else if (type_of_input == 1) {
		RandmOnCircle(p, n);
	}
	//for (int i = 0; i < n; i++) cout << p[i].x << "," << p[i].y << endl;

	while (collinear(p[0], p[1], p[2])) {
		int x = rand() % n;
		point temp = p[2];
		p[2] = p[x];
		p[x] = temp;
	}
	set<point> visiblep, visiblep1;
	facet t1,t2,t3,lt1,lt2;
	timer timer1,timer2;
	//确定每条线中的点顺时针排列
	//初始C,H,IC
	if (start(p[0], p[1], p[2])) {
		//cout << "check" << endl;
		t1 = { p[0],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t1)) {
				visiblep.insert(p[j]);
				//InversemapC.find(p[j])->second.insert(t);
			}
		}
		mapC.insert(pair<facet, set<point>>(t1, visiblep));
		H.insert(t1);
		visiblep.clear();
		t2 = { p[2],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t2)) {
				visiblep.insert(p[j]);
				//InversemapC.find(p[j])->second.insert(t);
			}
		}
		mapC.insert(pair<facet, set<point>>(t2, visiblep));
		H.insert(t2);
		visiblep.clear();
		t3 = { p[1],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t3)) {
				visiblep.insert(p[j]);
				//InversemapC.find(p[j])->second.insert(t);
			}
		}
		mapC.insert(pair<facet, set<point>>(t3, visiblep));
		H.insert(t3);
		visiblep.clear();
	}
	else {
		t1 = { p[0],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t1)) {
				visiblep.insert(p[j]);
				//InversemapC.find(p[j])->second.insert(t);
			}
		}
		mapC.insert(pair<facet, set<point>>(t1, visiblep));
		H.insert(t1);
		visiblep.clear();
		/*for (set<point>::iterator it = mapC[t].begin(); it != mapC[t].end(); it++)
		{
			cout << (*it).pivot << " 1 " << endl;
		}*/
		t2 = { p[1],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t2)) {
				visiblep.insert(p[j]);
				//InversemapC.find(p[j])->second.insert(t);
			}
		}
		mapC.insert(pair<facet, set<point>>(t2, visiblep));
		H.insert(t2);
		visiblep.clear();
		/*for (set<point>::iterator it = mapC[t].begin(); it != mapC[t].end(); it++)
		{
			cout << (*it).pivot << " 2 " << endl;
		}*/
		t3 = { p[2],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t3)) {
				visiblep.insert(p[j]);
				//InversemapC.find(p[j])->second.insert(t);
			}
		}
		mapC.insert(pair<facet, set<point>>(t3, visiblep));
		H.insert(t3);
		visiblep.clear();
	}
	set<facet> f;
	for (int j = 3; j < n; j++) {
		if (visible(p[j], t1)) {
			f.insert(t1);	
		}
		if (visible(p[j], t2)) {
			f.insert(t2);
		}
		if (visible(p[j], t3)) {
			f.insert(t3);
		}
		InversemapC.insert(pair<point, set<facet>>(p[j], f));
		f.clear();
	}

	//初始M
	facet test;
	facet test1;
	for (int i = 0; i < 3; i++) {
		//cout << i << endl;
		//cout << p[i].x << "," << p[i].y << endl;

		if (start(p[0], p[1], p[2])) {
			test = { p[i] ,p[(i + 2) % 3] };
			M.insert(pair<point, facet>(p[i], test));
			//cout << << endl;
			test = { p[(i + 1) % 3],p[i] };
			M.insert(pair<point, facet>(p[i], test));
			//cout <<  << endl;
		}
		else {
			test = { p[(i + 2) % 3], p[i] };
			M.insert(pair<point, facet>(p[i], test));
			//cout <<  << endl;
			test = { p[i], p[(i + 1) % 3] };
			M.insert(pair<point, facet>(p[i], test));
			//cout <<  << endl;
		}
	}

	/*int i=0;
	for (set<facet>::iterator it = InversemapC[p[5]].begin(); it != InversemapC[p[5]].end(); it++){
		cout << it->v1.pivot << " occurs " << endl;
		cout << i++ << endl;
	}*/
	point lr, rr;
	facet lt,rt,t;
	multimap<point, facet>::iterator miter;
	map<facet, set<point>>::iterator t1iter;
	map<facet, set<point>>::iterator t2iter;
	map<facet, set<point>>::iterator t3iter;
	set<point>::iterator piter;
	timer1.start();
	for (int i = 3; i < n; i++) {
		//cout << "-" << i << endl;
		

		//if(i == 16)for (set<facet>::iterator Hiter = H.begin(); Hiter != H.end(); ++Hiter) { cout << Hiter->v1.x << "," << Hiter->v1.y << " " << Hiter->v2.x << "," << Hiter->v2.y << endl; }
		//if (i == 4) { cout << InversemapC.find(p[i])->second.begin()->v1.pivot << endl; }
		//cout << InversemapC.size() << endl;
		if (InversemapC.find(p[i])->second.empty() || InversemapC.find(p[i]) == InversemapC.end()) {
			//cout << "countinue!" << endl;
			continue; 
		}
		set<facet> R = InversemapC[p[i]];
		//if (i == 16)cout << endl;
		//for (set<facet>::iterator Hiter = R.begin(); Hiter != R.end(); ++Hiter)
			//cout << Hiter->v1.x << "," << Hiter->v1.y << " " << Hiter->v2.x << "," << Hiter->v2.y << endl;
		
		
		
		//find r on boundary
		//lift
		for (set<facet>::iterator titer1 = R.begin(); titer1 != R.end(); ++titer1) {
			for (set<facet>::iterator titer2 = R.begin(); titer2 != R.end(); ) {
				if (titer1->v1.pivot == titer2->v2.pivot)break;
				if (++titer2 == R.end()) {
					lr = titer1->v1;
					lt2 = *titer1;
				}
			}
		}//find lr and t2
		for (set<facet>::iterator titer1 = R.begin(); titer1 != R.end(); ++titer1) {
			for (set<facet>::iterator titer2 = R.begin(); titer2 != R.end(); ) {
				if (titer1->v2.pivot == titer2->v1.pivot)break;
				if (++titer2 == R.end()) {
					rr = titer1->v2;
					t1 = *titer1;
				}
			}
		}
		if (rr.pivot == lt2.v2.pivot) {
			miter = M.find(lr);
			if (miter->second.v1.pivot == lt2.v1.pivot)lt1 = (++miter)->second;
			else lt1 = miter->second;
			miter = M.find(rr);
			if (miter->second.v1.pivot == t1.v1.pivot)t2 = (++miter)->second;
			else t2 = miter->second;
			lt = { lr,p[i] };
			t = { p[i], rr };
			//cout << "check1" << endl;
			t1iter = mapC.find(lt1);
			t2iter = mapC.find(lt2);
			t3iter = mapC.find(t2);
			piter = t2iter->second.begin();


			//C,IC insert1**************************
			for (int j = 0; j < t2iter->second.size(); j++, piter++) {
				InversemapC.find(*piter)->second.erase(lt2);
				if (visible((*piter), lt)) {
					visiblep.insert((*piter));
					InversemapC.find(*piter)->second.insert(lt);
				}
				if (visible((*piter), t)) {
					visiblep1.insert((*piter));
					InversemapC.find(*piter)->second.insert(t);
				}
			}
			//cout << "check2" << endl;
			piter = t1iter->second.begin();
			for (int j = 0; j < t1iter->second.size(); j++, piter++) {
				if (visible((*piter), lt)) {
					visiblep.insert((*piter));
					InversemapC.find(*piter)->second.insert(lt);
				}
			}
			mapC.erase(lt2);
			mapC.insert(pair<facet, set<point>>(lt, visiblep));
			
			//C,IC insert1*************************
			piter = t3iter->second.begin();
			for (int j = 0; j < t3iter->second.size(); j++, piter++) {
				if (visible((*piter), t)) {
					visiblep1.insert((*piter));
					InversemapC.find(*piter)->second.insert(t);
				}
			}
			mapC.insert(pair<facet, set<point>>(t, visiblep1));

			InversemapC.erase(InversemapC.find(p[i]));

			//cout << "t" << t2.v1.x << "," << t2.v1.y << " " << t2.v2.x << "," << t2.v2.y << endl;
			//cout << "t" << H.begin()->v1.x << "," << H.begin()->v1.y << " " << H.begin()->v2.x << "," << H.begin()->v2.y << endl;


			//H insert1,erase1**********************
			H.erase(lt2);
			//if(H.find(t2) == H.end())cout << "error" << endl;
			H.insert(lt);
			H.insert(t);
			//H insert1,erase1********************

			
			//M insert1,erase1*************************
			//if (i == 4)for (set<facet>::iterator Hiter = H.begin(); Hiter != H.end(); ++Hiter) cout << Hiter->v1.x << "," << Hiter->v1.y << " " << Hiter->v2.x << "," << Hiter->v2.y << endl;
			M.erase(lr);
			M.erase(rr);
			M.insert(pair<point, facet>(lr, lt1));
			M.insert(pair<point, facet>(lr, lt));
			M.insert(pair<point, facet>(lt.v2, lt));
			M.insert(pair<point, facet>(rr, t2));
			M.insert(pair<point, facet>(rr, t));
			M.insert(pair<point, facet>(t.v1, t));
			//cout << "cont" << endl;
			continue;
		}
		/*if (i == 4) {
			cout << "-----------" << endl;
			cout << lr.x << "," << lr.y << endl;
			cout << t2.v1.x << "," << t2.v1.y << " " << t2.v2.x << "," << t2.v2.y << endl;
		}*/

		//cout << "1" << endl;
		
		miter = M.find(lr);
		if (miter->second.v1.pivot == lt2.v1.pivot)lt1 = (++miter)->second;
		else lt1 = miter->second;
		lt = {lr,p[i]};
		//cout << "check1" << endl;
		t1iter = mapC.find(lt1);
		t2iter = mapC.find(lt2);
		piter = t2iter->second.begin();


		//C,IC insert1**************************
		for (int j = 0; j < t2iter->second.size(); j++, piter++) {
			InversemapC.find(*piter)->second.erase(lt2);
			if (visible((*piter), lt)) {
				visiblep.insert((*piter));
				InversemapC.find(*piter)->second.insert(lt);
			}
		}
		//cout << "check2" << endl;
		piter = t1iter->second.begin();
		for (int j = 0; j < t1iter->second.size(); j++, piter++) {
			if (visible((*piter), lt)) {
				visiblep.insert((*piter));
				InversemapC.find(*piter)->second.insert(lt);
			}
		}
		//cout << "2" << endl;
		mapC.erase(lt2);
		mapC.insert(pair<facet, set<point>>(lt, visiblep));
		//C,IC insert1*************************




		//cout << "t" << t2.v1.x << "," << t2.v1.y << " " << t2.v2.x << "," << t2.v2.y << endl;
		//cout << "t" << H.begin()->v1.x << "," << H.begin()->v1.y << " " << H.begin()->v2.x << "," << H.begin()->v2.y << endl;
		
		
		//H insert1,erase1**********************
		H.erase(lt2);
		//if(H.find(t2) == H.end())cout << "error" << endl;
		H.insert(lt);
		//H insert1,erase1********************


		//M insert1,erase1*************************
		//if (i == 4)for (set<facet>::iterator Hiter = H.begin(); Hiter != H.end(); ++Hiter) cout << Hiter->v1.x << "," << Hiter->v1.y << " " << Hiter->v2.x << "," << Hiter->v2.y << endl;
		M.erase(lr);
		M.insert(pair<point, facet>(lr, lt1));
		M.insert(pair<point, facet>(lr, lt));
		M.insert(pair<point, facet>(lt.v2, lt));
		//M insert1,erase1***********************


		//right
		/*for (set<facet>::iterator titer1 = R.begin(); titer1 != R.end(); ++titer1) {
			for (set<facet>::iterator titer2 = R.begin(); titer2 != R.end(); ) {
				if (titer1->v2.pivot == titer2->v1.pivot)break;
				if (++titer2 == R.end()) {
					rr = titer1->v2;
					t1 = *titer1;
				}
			}
		}
		
		for (multimap<point, facet>::iterator m = M.begin(); m != M.end();++m) {
				cout << m->first.pivot << "--" << m->second.v1.pivot << "," << m->second.v2.pivot << endl;
		}*/
		
		
		miter = M.find(rr);
		if (miter->second.v1.pivot == t1.v1.pivot)t2 = (++miter)->second;
		else t2 = miter->second;
		t = { p[i], rr};

		t1iter = mapC.find(t1);
		t2iter = mapC.find(t2);
		piter = t1iter->second.begin();

		//C,IC insert2*************************
		for (int j = 0; j < t1iter->second.size(); j++, piter++) {
			InversemapC.find(*piter)->second.erase(t1);
			if (visible((*piter), t)) {
				visiblep.insert((*piter));
				InversemapC.find(*piter)->second.insert(t);
			}
		}//cout << "check1" << endl;
		//cout << rr.pivot << endl;
		for (piter = t2iter->second.begin(); piter != t2iter->second.end(); ++piter) {
			if (visible((*piter), t)) {
				visiblep.insert((*piter));
				InversemapC.find(*piter)->second.insert(t);
			}
		}//cout << "check2" << endl;
		
		mapC.erase(t1);
		mapC.insert(pair<facet, set<point>>(t, visiblep));
		//C,IC insert2*************************
		//cout << "check3" << endl;
		//H insert2**********************
		H.erase(t1);
		H.insert(t);
		//H insert2**********************
		M.erase(rr);
		M.insert(pair<point, facet>(rr, t2));
		M.insert(pair<point, facet>(rr, t));
		M.insert(pair<point, facet>(t.v1, t));
		//cout << "3" << endl;
		//IC erase1*************************
		InversemapC.erase(InversemapC.find(p[i]));
		//IC erase1*************************
		//cout << "check4" << endl;
		R.erase(t1);
		R.erase(lt2);
		if (R.empty())continue;
		set<facet>::iterator titer = R.begin();
		for (int j = 0; j < R.size(); j++, titer++) {
			map<point, set<facet>>::iterator ICiter = InversemapC.begin();
			for (int k = 0; k < InversemapC.size(); k++, ICiter++) {
				ICiter->second.erase(*titer);
			}
			H.erase(*titer);
			mapC.erase(*titer);
			M.erase(titer->v2);
		}//cout << "check5" << endl;
		
	}
	timer1.stop();
	
	timer2.start();
	memset(s, 0, sizeof(s));
	sort(p, p + n, cmp1);//排序找出纵坐标最小的值 
	//for (int i = 0; i < n; i++) cout << p[i].x << "," << p[i].y << endl;
	s[0] = p[0];
	sort(p + 1, p + n, cmp2);//极角排序
	s[1] = p[1];//找到p1点 
	int top = 1;
	for (int i = 2; i < n; i++) {
		while (cross_product(s[top - 1], s[top], p[i]) <= 0) {
			top--;//向右转，这个中间点删除
		}
		s[++top] = p[i];//向左转，添加 
	}
	timer2.stop();
	sort(s, s + top + 1, cmp1);
	//		for (int i = 0; i < top+1; i++) cout << s[i].x << "," << s[i].y << endl;
	for (int i = 0; i < top; i++) {
		//			cout << i << ":" << top << endl;
		while (s[i].x == s[i + 1].x && s[i].y == s[i + 1].y) {
			for (int j = i + 1; j < top; j++) {
				s[j] = s[j + 1];
			}
			top--;
		}
	}
	if (H.size() >= top + 1) {
		for (int i = 0; i < top + 1; i++) {
			set<facet>::iterator check = H.begin();
			for (int j = 0; j < H.size(); j++, check++) {
				if (s[i].pivot == check->v1.pivot) break;
				if (j == H.size() - 1)return 0;
			}
		}
		cout << "correct" << endl;
	}/**/
	cout << "time1: " << timer1.get_total() << endl;
	cout << "time2: " << timer2.get_total() << endl;
	//for (set<facet>::iterator Hiter = H.begin(); Hiter != H.end(); ++Hiter) { cout << Hiter->v1.x << "," << Hiter->v1.y << " " << Hiter->v2.x << "," << Hiter->v2.y << endl; }

	return 0;


}