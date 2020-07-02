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

const  double PI = 3.1415926;
const  double EPSINON = 0.0000001;
int count1 = 0;
int count2 = 0;
int count3 = 0;
int count4 = 0;
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
	int pivot;
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

namespace std {
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
}
using namespace std;

unordered_map<facet, vector<point>> mapC;
unordered_map<point, facet> M;
unordered_set<Facet> H;


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

Facet GetValue(point r, facet t) {

}

void ProcessRidge(facet t1, point r, facet t2) {
	//testtime1.start();
	unordered_map<facet, vector<point>>::iterator t1iter;
	unordered_map<facet, vector<point>>::iterator t2iter;
	vector<point> vp;
	
	//cout << "check1" << endl;
	t1iter = mapC.find(t1);
	t2iter = mapC.find(t2);
	//testtime1.stop();
	if (t1iter->second.empty() && t2iter->second.empty()) {
		return;
	}
	else if (!(t2iter->second.empty()) &&  !(t1iter->second.empty()) && t1iter->second.back() == t2iter->second.back()) {
		H.erase(t1);
		H.erase(t2);
	}
	else if (t1iter->second.empty() || (t2iter->second.size() != 0 && (t2iter->second.back() < t1iter->second.back()))) {//左
		//testtime2.start();
		point p = (t2iter->second.back());
		point temp;
		facet t;
		t.v1 = r;
		t.v2 = p;
		//cout << "3check2" << endl;
		cilk_for (int i = 0; i < t2iter->second.size(); i++) {
			//cout << i << endl;
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
		//cout << "3check3" << endl;
		cilk_for (int i = 0; i < t1iter->second.size(); i++) {
			//cout << i << endl;
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
		}
		
		//testtime2.stop();
		//testtime3.start();
		mapC.insert(pair<facet, vector<point>>(t, vp));

		
		H.erase(t2);
		H.insert(t);
		point tempr;
		facet tempt; 
		//testtime3.stop();
		for (int i = 0; i < 2; i++) {
			
			if (i) tempr = t.v1;
			else tempr = t.v2;
			if (r.pivot == tempr.pivot) {
				ProcessRidge(t1, r, t);
			}
			else{
				//testtime4.start();
				unordered_map<point, facet>::iterator it = M.find(tempr);
				//testtime4.stop();
				if (it == M.end()) {
					//testtime4.start();
					M.insert(pair<point, facet>(tempr, t));
					//testtime4.stop();
				}
				else {
					ProcessRidge(t, tempr, it->second);
				}
			}
		}
	}
	else {//右
		//testtime5.start();
		point p = t1iter->second.back();
		point temp;
		facet t;
		t.v1 = p;
		t.v2 = r;
		//cout << "4check2" << endl;
		for (int i = 0; i < t2iter->second.size(); i++) {
			//cout << i << endl;
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
		//cout << "4check3" << endl;
		for (int i = 0; i < t1iter->second.size(); i++) {
			//cout << i << endl;
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
		}
		//cout << "4check4" << endl;
		mapC.insert(pair<facet, vector<point>>(t, vp));


		H.erase(t1);
		H.insert(t);
		point tempr;
		facet tempt;
		//testtime5.stop();
		for (int i = 0; i < 2; i++) {
			//testtime5.start();
			if (i) tempr = t.v1;
			else tempr = t.v2;
			//testtime5.stop();
			if (r.pivot == tempr.pivot) {
				ProcessRidge(t, r, t2);
			}
			else {
				//testtime5.start();
				unordered_map<point, facet>::iterator it = M.find(tempr);
				//testtime5.stop();
				if (it == M.end()) {
					//testtime5.start();
					M.insert(pair<point, facet>(tempr, t));
					//testtime5.stop();
				}
				else {
					ProcessRidge(it->second, tempr, t);
				}
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

	while (collinear(p[0], p[1], p[2])) {//排除初始三点共线
		int x = rand() % n;
		point temp = p[2];
		p[2] = p[x];
		p[x] = temp;
	}
	vector<point> visiblep;
	point temp;
	facet t, facet1, facet2;
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
		sort(visiblep.begin(), visiblep.end(), sortfunction);
		mapC.insert(pair<facet, vector<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();
		t = { p[2],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.push_back(p[j]);
			}
		}
		sort(visiblep.begin(), visiblep.end(), sortfunction);
		mapC.insert(pair<facet, vector<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();
		t = { p[1],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.push_back(p[j]);
			}
		}
		sort(visiblep.begin(), visiblep.end(), sortfunction);
		mapC.insert(pair<facet, vector<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();
	}
	else {
		t = { p[0],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.push_back(p[j]);
			}
		}
		sort(visiblep.begin(), visiblep.end(), sortfunction);
		mapC.insert(pair<facet, vector<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();

		t = { p[1],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.push_back(p[j]);
			}
		}
		sort(visiblep.begin(), visiblep.end(), sortfunction);
		mapC.insert(pair<facet, vector<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();

		t = { p[2],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.push_back(p[j]);
			}
		}
		sort(visiblep.begin(), visiblep.end(), sortfunction);
		mapC.insert(pair<facet, vector<point>>(t, visiblep));
		H.insert(t);
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
	//cout << "find times " << count1 << endl;
	/*cout << endl;
	cout << "init and find t1,t2:			" << testtime1.get_total() << endl;
	//cout << "find times " << count2 << endl;
	
	cout << "check visible:				" << testtime2.get_total() << endl;
	//cout << "find times " << count3 << endl;
	cout << "mapC insert and H insert and erase:	" << testtime3.get_total() << endl;
	//cout << "find times " << count4 << endl;
	cout << "M map check:				" << testtime4.get_total() << endl;
	cout << "case 4 total time:			" << testtime5.get_total() << endl;
	*/
	return 0;
}
