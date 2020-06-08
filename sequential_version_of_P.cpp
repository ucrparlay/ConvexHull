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
int count1 = 0;


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
}s[10000];//

typedef struct facet
{
	point v1, v2;
	bool operator < (const facet &o) const
	{
		if (v1.pivot == o.v1.pivot)
			return v2.pivot < o.v2.pivot;
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
multimap<point, facet> M;
set<Facet> H;
//hash_map<int, ridgefacet> R;

void RandmOnCircle(point *p, int n) {
	for (int i = 0; i < n; i++) {
		p[i].x = (double)n + n*cos((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
		p[i].y = (double)n + n*sin((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
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
	if ((v1.x - v3.x) * (v2.y - v3.y) - (v1.y - v3.y) * (v2.x - v3.x) < EPSINON  && (v1.x - v3.x) * (v2.y - v3.y) - (v1.y - v3.y) * (v2.x - v3.x) > -EPSINON) return true;
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
			return (atan2(b.y - s[0].y, b.x - s[0].x)-(atan2(a.y - s[0].y, a.x - s[0].x)) > EPSINON);
		return (b.x - a.x > EPSINON) ;
}

bool mins(point a, point b) {
	return a.pivot < b.pivot;
}

/*bool InsertandSet(point r, facet t) {
	multimap<point, facet>::iterator it = M.find(r);
	if (it->second == NULL) M.insert(pair<point, facet>(r, t));
	else {
		it->second
	}
}*/

Facet GetValue(point r, facet t) {

}/**/

void ProcessRidge(facet t1, point r, facet t2) {
	map<facet, set<point>>::iterator t1iter;
	map<facet, set<point>>::iterator t2iter;
	set<point>::iterator piter;
	set<point> visiblep;
	//if (++count1 >= 25)return;
	/*for (map<facet, set<point>>::iterator Hiter = mapC.begin(); Hiter != mapC.end(); ++Hiter) {
		cout << "facet:" << Hiter->first.v1.x << "," << Hiter->first.v1.y << " " << Hiter->first.v2.x << "," << Hiter->first.v2.y << " " << Hiter->first.v2.pivot << endl;
		for (set<point>::iterator Hiter1 = Hiter->second.begin(); Hiter1 != Hiter->second.end(); ++Hiter1) { cout << Hiter1->x << "," << Hiter1->y << " " << Hiter1->pivot << endl; }

	}*/

	t1iter = mapC.find(t1);
	//for (set<point>::iterator it = t1iter->second.begin(); it != t1iter->second.end(); ++it)cout << "t1 visible:" << it->pivot << endl;
	t2iter = mapC.find(t2);
	//for (set<point>::iterator it = t2iter->second.begin(); it != t2iter->second.end(); ++it)cout << "t2 visible:" << it->pivot << endl;
	//cout << "check ProcessRidge" << endl;
	/*if (t1iter == mapC.end()) {
		cout << "1" << endl;
	}
	if (t2iter == mapC.end()) {
		cout << "2" << endl;
	}*/
	
	if (t1iter->second.size() == 0 && t2iter->second.size() == 0) {
		//cout << "check ProcessRidge 1" << endl;
		return;
	}
	else if (t1iter->second.begin()->pivot == t2iter->second.begin()->pivot) {
		//cout << "check ProcessRidge 2" << endl;
		H.erase(t1);
		H.erase(t2);
		//mapC.erase(t1iter);
		//mapC.erase(t2iter);
		//cout << t1iter->second.begin()->pivot << " " << t2iter->second.begin()->pivot  << endl;
	}
	else if (((t1iter->second.begin()->pivot > t2iter->second.begin()->pivot) && t2iter->second.size() != 0) || t1iter->second.size() == 0) {//左
		//cout << "check ProcessRidge 3" << endl;
		//cout << "M : " << M.size() << endl;
		point p = (*t2iter->second.begin());
		facet t;
		t.v1 = r;
		t.v2 = p;
		piter = t2iter->second.begin();
		//cout << (*piter).pivot << " " << t2iter->second.size() << endl;
		//C(t)***************************************
		if(t1iter == mapC.end() || t1iter->second.size() == 0){
			for (int i = 0; i < t2iter->second.size(); i++, piter++) {
				if (visible((*piter), t))visiblep.insert((*piter));
			}
		}
		else {
			for (int i = 0; i < t2iter->second.size(); i++, piter++) {
				//if (piter->pivot == t1iter->second.begin()->pivot)break;
				if (visible((*piter), t))visiblep.insert((*piter));
			}
			for (piter = t1iter->second.begin(); piter != t1iter->second.end(); ++piter) {
				if (visible((*piter), t))visiblep.insert((*piter));
			}
		}
		
		//mapC.erase(t2iter);


		/*for (set<point>::iterator it = visiblep.begin(); it != visiblep.end(); ++it)
			cout << "t1 visible:" << it->pivot << endl;
		cout << endl;
		for (map<facet, set<point>>::iterator Hiter = mapC.begin(); Hiter != mapC.end(); ++Hiter) {
			cout << "facet:" << Hiter->first.v1.x << "," << Hiter->first.v1.y << " " << Hiter->first.v2.x << "," << Hiter->first.v2.y << " " << Hiter->first.v2.pivot << endl;
			for (set<point>::iterator Hiter1 = Hiter->second.begin(); Hiter1 != Hiter->second.end(); ++Hiter1) { cout << Hiter1->x << "," << Hiter1->y << " " << Hiter1->pivot << endl; }

		}*/
		mapC.insert(pair<facet, set<point>>(t, visiblep));
		visiblep.clear();
		/*for (map<facet, set<point>>::iterator Hiter = mapC.begin(); Hiter != mapC.end(); ++Hiter) {
			cout << "facet:" << Hiter->first.v1.x << "," << Hiter->first.v1.y << " " << Hiter->first.v2.x << "," << Hiter->first.v2.y << " " << Hiter->first.v2.pivot << endl;
			for (set<point>::iterator Hiter1 = Hiter->second.begin(); Hiter1 != Hiter->second.end(); ++Hiter1) { cout << Hiter1->x << "," << Hiter1->y << " " << Hiter1->pivot << endl; }

		}
		for (set<point>::iterator it = mapC.find(t)->second.begin(); it != mapC.find(t)->second.end(); ++it)
			cout << "t1 visible:" << it->pivot << endl;
		cout << endl;*/


		//H**********************************************
		H.erase(t2);
		H.insert(t);
		point tempr;
		facet tempt;
		//r,r'*********************************************
		for (int i = 0; i < 2; i++) {
			if (i) tempr = t.v1;
			else tempr = t.v2;
			if (r.pivot == tempr.pivot) {
				//cout << "3 R" << endl;
				//cout << r.pivot << endl;
				//cout << t.v1.pivot << "," << t.v2.pivot << endl;
				ProcessRidge(t1, r, t);
			}
			else{
				//for (multimap<point, facet>::iterator it = M.begin(); it != M.end(); ++it)cout << "M:" << it->first.pivot << " " << it->second.v1.pivot << "," << it->second.v2.pivot << endl;
				//cout << endl;
				multimap<point, facet>::iterator it = M.find(tempr);
				//cout << tempr.pivot << endl;
				//cout << it->second.v1.pivot << "," << it->second.v2.pivot << endl;
				if (it == M.end()) {
					//cout << "3 if" << endl;
					M.insert(pair<point, facet>(tempr, t));
				}
				else {
					//cout << "3 else" << endl;
					ProcessRidge(t, tempr, it->second);
				}
			}
		}
	}
	else {//右
		//cout << "check ProcessRidge 4" << endl;
		point p = (*t1iter->second.begin());
		facet t;
		t.v1 = p;
		t.v2 = r;
		piter = t1iter->second.begin();
		//cout << (*piter).pivot << " " << t2iter->second.size() << endl;
		if (t2iter == mapC.end() || t2iter->second.size() == 0) {
			//cout << "check" << endl;
			for (int i = 0; i < t1iter->second.size(); i++, piter++) {
				if (visible((*piter), t))visiblep.insert((*piter));
				//if (t1iter->second.begin()->pivot == 7) cout << i << endl;
			}
		}
		else {
			for (int i = 0; i < t1iter->second.size(); i++, piter++) {
				//if (piter->pivot == t2iter->second.begin()->pivot)break;
				if (visible((*piter), t))visiblep.insert((*piter));
			}
			for (piter = t2iter->second.begin(); piter != t2iter->second.end(); piter++) {
				if (visible((*piter), t))visiblep.insert((*piter));
			}
		}
		//mapC.erase(t1iter);
		/*if (t1iter->second.begin()->pivot == 7) {
			for (set<point>::iterator it = visiblep.begin(); it != visiblep.end(); ++it)
				cout << "t visible:" << it->pivot << endl;
			cout << endl;
		}*/
		
		mapC.insert(pair<facet, set<point>>(t, visiblep));
		visiblep.clear();

		H.erase(t1);
		H.insert(t);
		point tempr;
		facet tempt;
		for (int i = 0; i < 2; i++) {
			if (i) tempr = t.v1;
			else tempr = t.v2;
			if (r.pivot == tempr.pivot) {
				//cout << "4 R" << endl;
				//cout << r.pivot << endl;
				//cout << t.v1.pivot << "," << t.v2.pivot << endl;
				ProcessRidge(t, r, t2);
			}
			else {
				multimap<point, facet>::iterator it = M.find(tempr);
				if (it == M.end()) {
					//cout << "4 if" << endl;
					M.insert(pair<point, facet>(tempr, t));
				}
				else {
					//cout << "4 else" << endl;
					ProcessRidge(it->second, tempr, t);
				}
			}
		}
	}
	//cout << "end" << endl;
	return;
}/**/

int main(int argc, char** argv) {
	
	if (argc != 3) {
		cout << "Usage: ./qsort [num_elements] [type of input]" << endl;
		return 0;
	}
	int n = atoi(argv[1]);
	int type_of_input = atoi(argv[2]);
	point* p = new point[n];
	
	if(type_of_input == 0){
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
	set<point> visiblep;
	facet t, facet1, facet2;
	//确定每条线中的点顺时针排列
	//初始C,H
	if (start(p[0], p[1], p[2])) {
		//cout << "check" << endl;
		t = { p[0],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, set<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();
		t = { p[2],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, set<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();
		t = { p[1],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, set<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();
	}
	else {
		t = { p[0],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, set<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();
		/*for (set<point>::iterator it = mapC[t].begin(); it != mapC[t].end(); it++)
		{
			cout << (*it).pivot << " 1 " << endl;
		}*/
		t = { p[1],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, set<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();
		/*for (set<point>::iterator it = mapC[t].begin(); it != mapC[t].end(); it++)
		{
			cout << (*it).pivot << " 2 " << endl;
		}*/
		t = { p[2],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, set<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();
	}
	/*for (set<point>::iterator it = mapC[t = { p[0],p[1] }].begin(); it != mapC[t = { p[0],p[1] }].end(); ++it)
	{
		cout << (*it).pivot << " occurs " << endl;
	}
	cout << endl;
	for (set<point>::iterator it = mapC[t = { p[1],p[2] }].begin(); it != mapC[t = { p[1],p[2] }].end(); ++it)
	{
		cout << (*it).pivot << " occurs " << endl;
	}
	cout << endl;
	for (set<point>::iterator it = mapC[t = { p[2],p[0] }].begin(); it != mapC[t = { p[2],p[0] }].end(); ++it)
	{
		cout << (*it).pivot << " occurs " << endl;
	}
	cout << endl;
	

	cout << "check CH" << endl;*/
	//timer t1,t2,t3,t4,t5;
	timer timer1, timer2, timer3, timer4, timer5, timer6, timer7;
	facet test;
	facet test1;
	//cout << sizeof(p) << endl;
	/*for (int i = 0; i < 3; i++) {
		set<point> visiblep;
		//mapC.insert(pair<int, vector<point>>(i, visiblep));
		
		for (int j = 3; j < n; j++) {
			if (visible(p[j], H[i])) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, set<point>>(H[i], visiblep));
	}

	for (int i = 0; i < mapC.size(); i++) 
		for (set<point>::iterator it = mapC.find(H[i]); it != numSet.end(); ++it)
		{
			cout << *it << " occurs " << endl;
		}*/
	//初始M
	bool use;
	int flag;

	for (int i = 0; i < 3; i++) {
		//cout << i << endl;
		//cout << p[i].x << "," << p[i].y << endl;
		if (start(p[0], p[1], p[2])) {
			test = {p[i] ,p[(i+2)%3] };
			M.insert(pair<point, facet>(p[i], test));
			//cout << << endl;
			test = { p[(i+1)%3],p[i] };
			M.insert(pair<point, facet>(p[i], test));
			//cout <<  << endl;
		}
		else {
			test = { p[(i + 2) % 3], p[i] };
			M.insert(pair<point, facet>(p[i], test));
			//cout <<  << endl;
			test = { p[i], p[(i + 1) % 3]};
			M.insert(pair<point, facet>(p[i], test));
			//cout <<  << endl;
		}
	}
	timer1.start();
	//cout << "check M" << endl;
	if (start(p[0], p[1], p[2])) {
		ProcessRidge(facet1 = { p[1],p[0] }, p[0], facet2 = { p[0],p[2] });
		ProcessRidge(facet1 = { p[2],p[1] }, p[1], facet2 = { p[1],p[0] });
		ProcessRidge(facet1 = { p[0],p[2] }, p[2], facet2 = { p[2],p[1] });
	}
	else {
		ProcessRidge(facet1 = { p[2],p[0] }, p[0], facet2 = { p[0],p[1] });
		//cout << "1" << endl;
		ProcessRidge(facet1 = { p[0],p[1] }, p[1], facet2 = { p[1],p[2] });
		//cout << "2" << endl;
		//
		ProcessRidge(facet1 = { p[1],p[2] }, p[2], facet2 = { p[2],p[0] });if (count1 >= 10)return 0;
	}
	timer1.stop();
	//cout << "check output" << endl;
	/*for (set<facet>::iterator Hiter = H.begin(); Hiter != H.end(); ++Hiter) {
		cout << Hiter->v1.x << "," << Hiter->v1.y << " " << Hiter->v2.x << "," << Hiter->v2.y << endl;
	}*/

	//初始R
	//R.insert(pair<int, ridgefacet>(hash1(p[0].pivot), test1));
	//M.insert(pair<point, facet>(p[1], test = { p[2],p[0] }));
	//for (int i = 0; i < mapC[1].size(); i++) cout << mapC[1][i].x << endl;
	//for (int i = 0; i < M.size(); i++) cout << M.find(H[i].v2). << endl;




	/*for (int i = 0; i < M.size(); i++) {
		cout << i << endl;
		multimap<point, facet>::iterator m;
		m = M.find(p[i]);
		if(m == M.end())cout << i << endl;
		for (int k = 0; k != M.count(p[i]); k++, m++)
			cout << m->first.x << "," << m->first.y << "--" << m->second.v1.x << "," << m->second.v1.y
			<< " " << m->second.v2.x << "," << m->second.v2.y << endl;
		cout << m->first.x << "," << m->first.y << endl;
		cout << m->second[0].v1.x << "," << m->second[0].v1.y
			<< " " << m->second[0].v2.x << "," << m->second[0].v2.y
			<< "--" << m->second[1].v1.x << "," << m->second[1].v1.y
			<< " " << m->second[1].v2.x << "," << m->second[1].v2.y << endl;
	}*/

	/*timer2.start();
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
	}*/
	/*for (set<facet>::iterator Hiter = H.begin(); Hiter != H.end(); ++Hiter) { cout << Hiter->v1.x << "," << Hiter->v1.y << " " << Hiter->v2.x << "," << Hiter->v2.y << endl; }
	cout << endl;
	for (int i = 0; i < top + 1; i++) cout << s[i].x << "," << s[i].y << endl;
	if (H.size() >= top + 1) {
		for (int i = 0; i < top + 1; i++) {
			set<facet>::iterator check = H.begin();
			for (int j = 0; j < H.size(); j++, check++) {
				if (s[i].pivot == check->v1.pivot) break;
				if (j == H.size() - 1)return 0;
			}
		}
		cout << "correct" << endl;
	}*/
//	for (int i = 0; i < H.size(); i++) cout << H[i].v1.x << " " << H[i].v1.y << " " << H[i].v2.x << " " << H[i].v2.y << endl;
//for (int i = 0; i < top+1; i++) cout << s[i].x << "," << s[i].y << endl;


		/*if (H.size() >= top + 1) {
			cout << "size1: " << H.size() << endl;
			cout << "size2: " << top + 1 << endl;
			cout << "time1: " << t1.get_total() << endl;
			cout << "time2: " << t2.get_total() << endl;
			cout << "time3: " << t3.get_total() << endl;
			//cout << "time4: " << t4.get_total() << endl;
			//cout << "time5: " << t5.get_total() << endl;
		}
		else {
			for (int i = 0; i < top + 1; i++) {
				for (int j = 0; j < H.size(); j++) {
					if (s[i].x == H[j].v1.x && s[i].y == H[j].v1.y) break;
					if (j == H.size() - 1)cout << s[i].x << "," << s[i].y << endl;;
				}
			}
			cout << "size1: " << H.size() << endl;
			cout << "size2: " << top + 1 << endl;
			cout << "time1: " << t1.get_total() << endl;
			cout << "time2: " << t2.get_total() << endl;
			cout << "time3: " << t3.get_total() << endl;
		}
	*/
	cout << "total time: " << timer1.get_total() << endl;
	//cout << "jarvis time: " << timer2.get_total() << endl;
	cout << "size1: " << H.size() << endl;
	//cout << "size2: " << top + 1 << endl;

	return 0;
}
