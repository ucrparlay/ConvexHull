#pragma once
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <atomic>
#include <mutex>
#include <string.h>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <set>
#include "dataset.h"
#include "CheckVisible.h"

const int MAX = 100000000;
namespace std {
	template<>
	struct hash<Facet> {
		std::size_t operator()(const Facet& o) const {
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

unordered_map<facet, unordered_set<point>> mapC;
unordered_map<point, facet> M;
unordered_set<Facet> H;
void init(point *p, int n){
	while (collinear(p[0], p[1], p[2])) {//排除初始三点共线
		int x = rand() % n;
		point temp = p[2];
		p[2] = p[x];
		p[x] = temp;
	}
	unordered_set<point> visiblep;
	facet t;
	//确定每条线中的点顺时针排列
	//初始C,H
	if (start(p[0], p[1], p[2])) {
		t = { p[0],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, unordered_set<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();
		t = { p[2],p[1] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, unordered_set<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();
		t = { p[1],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, unordered_set<point>>(t, visiblep));
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
		mapC.insert(pair<facet, unordered_set<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();

		t = { p[1],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, unordered_set<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();

		t = { p[2],p[0] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, unordered_set<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();
	}
	facet tempf;
	//初始M
	/*for (int i = 0; i < 3; i++) {

		if (start(p[0], p[1], p[2])) {
			tempf = { p[i] ,p[(i + 2) % 3] };
			M.insert(pair<point, facet>(p[i], tempf));

			tempf = { p[(i + 1) % 3],p[i] };
			M.insert(pair<point, facet>(p[i], tempf));
		}
		else {
			tempf = { p[(i + 2) % 3], p[i] };
			M.insert(pair<point, facet>(p[i], tempf));

			tempf = { p[i], p[(i + 1) % 3] };
			M.insert(pair<point, facet>(p[i], tempf));
		}
	}*/
}

Facet GetValue(point r, facet t) {

}

void ProcessRidge(facet t1, point r, facet t2) {
	//testtime3.start();
	unordered_map<facet, set<point>>::iterator t1iter;
	unordered_map<facet, set<point>>::iterator t2iter;
	unordered_set<point>::iterator piter;
	unordered_set<point> visiblep;
	point t1min,t2min; 
	
	t1iter = mapC.find(t1);
	t2iter = mapC.find(t2);
	//testtime3.stop();
	if (t1iter->second.size() == 0 && t2iter->second.size() == 0) {
		return;
	}
	else {
		//testtime1.start();
		t2min = {0,0,MAX};
		t1min = { 0,0,MAX };
		for (unordered_set<point>::iterator it = t1iter->second.begin(); it != t1iter->second.end();++it) {
			if (*it < t1min)t1min = *it;
		}
		for (unordered_set<point>::iterator it = t2iter->second.begin(); it != t2iter->second.end(); ++it) {
			if (*it < t2min) {
				t2min = *it;
			}
		}//testtime1.stop();
		if (t2min.pivot == t1min.pivot) {
		H.erase(t1);
		H.erase(t2);
		}
		else if (((t2min.pivot < t1min.pivot) && t2iter->second.size() != 0) || t1iter->second.size() == 0) {//左
			//testtime2.start();
			point p = t2min;
			facet t;
			t.v1 = r;
			t.v2 = p;
			piter = t2iter->second.begin();
			if (t1iter == mapC.end() || t1iter->second.size() == 0) {
				for (int i = 0; i < t2iter->second.size(); i++, piter++) {
					if (visible((*piter), t))visiblep.insert((*piter));
				}
			}
			else {
				for (int i = 0; i < t2iter->second.size(); i++, piter++) {
					if (visible((*piter), t))visiblep.insert((*piter));
				}
				for (piter = t1iter->second.begin(); piter != t1iter->second.end(); ++piter) {
					if (visible((*piter), t))visiblep.insert((*piter));
				}
			}
			//testtime2.stop();
			//testtime5.start();
			mapC.insert(pair<facet, unordered_set<point>>(t, visiblep));

			H.erase(t2);
			H.insert(t);
			point tempr;
			facet tempt;
			//testtime5.stop();
			for (int i = 0; i < 2; i++) {
				if (i) tempr = t.v1;
				else tempr = t.v2;
				if (r.pivot == tempr.pivot) {
					ProcessRidge(t1, r, t);
				}
				else {
					//testtime6.start();
					unordered_map<point, facet>::iterator it = M.find(tempr);
					//testtime6.stop();
					if (it == M.end()) {
						//testtime6.start();
						M.insert(pair<point, facet>(tempr, t));
						//testtime6.stop();
					}
					else {
						ProcessRidge(t, tempr, it->second);
					}
				}
			}
		}
		else {//右
			//testtime4.start();
			point p = t1min;
			facet t;
			t.v1 = p;
			t.v2 = r;
			piter = t1iter->second.begin();
			if (t2iter == mapC.end() || t2iter->second.size() == 0) {
				for (int i = 0; i < t1iter->second.size(); i++, piter++) {
					if (visible((*piter), t))visiblep.insert((*piter));
				}
			}
			else {
				for (int i = 0; i < t1iter->second.size(); i++, piter++) {
					if (visible((*piter), t))visiblep.insert((*piter));
				}
				for (piter = t2iter->second.begin(); piter != t2iter->second.end(); piter++) {
					if (visible((*piter), t))visiblep.insert((*piter));
				}
			}

			mapC.insert(pair<facet, unordered_set<point>>(t, visiblep));
			
			H.erase(t1);
			H.insert(t);
			point tempr;
			facet tempt;
			//testtime4.stop();
			for (int i = 0; i < 2; i++) {
				if (i) tempr = t.v1;
				else tempr = t.v2;
				if (r.pivot == tempr.pivot) {
					ProcessRidge(t, r, t2);
				}
				else {
					//testtime4.start();
					unordered_map<point, facet>::iterator it = M.find(tempr);
					//testtime4.stop();
					if (it == M.end()) {
						//testtime4.start();
						M.insert(pair<point, facet>(tempr, t));
						//testtime4.stop();
					}
					else {
						ProcessRidge(it->second, tempr, t);
					}
				}
			}
		}
	}
	return;
}
