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
namespace std {
	template<>
	struct hash<Facet> {
		std::size_t operator()(const Facet& o) const {
			using std::hash;
			return (hash<int>{}(o.v1.pivot) ^ (hash<int>{}(o.v2.pivot) << 1));
		}
	};
}
using namespace std;

unordered_map<facet, set<point>> mapC;
multimap<point, facet> M;
set<Facet> H;
void init(point *p, int n){
	while (collinear(p[0], p[1], p[2])) {//排除初始三点共线
		int x = rand() % n;
		point temp = p[2];
		p[2] = p[x];
		p[x] = temp;
	}
	set<point> visiblep;
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

		t = { p[1],p[2] };
		for (int j = 3; j < n; j++) {
			if (visible(p[j], t)) {
				visiblep.insert(p[j]);
			}
		}
		mapC.insert(pair<facet, set<point>>(t, visiblep));
		H.insert(t);
		visiblep.clear();

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
	facet tempf;
	//初始M
	for (int i = 0; i < 3; i++) {

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
	}
}

Facet GetValue(point r, facet t) {

}

void ProcessRidge(facet t1, point r, facet t2) {
	unordered_map<facet, set<point>>::iterator t1iter;
	unordered_map<facet, set<point>>::iterator t2iter;
	set<point>::iterator piter;
	set<point> visiblep;
	t1iter = mapC.find(t1);
	t2iter = mapC.find(t2);

	if (t1iter->second.size() == 0 && t2iter->second.size() == 0) {
		return;
	}
	else if (t1iter->second.begin()->pivot == t2iter->second.begin()->pivot) {
		H.erase(t1);
		H.erase(t2);
	}
	else if (((t1iter->second.begin()->pivot > t2iter->second.begin()->pivot) && t2iter->second.size() != 0) || t1iter->second.size() == 0) {//左
		point p = (*t2iter->second.begin());
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

		mapC.insert(pair<facet, set<point>>(t, visiblep));
		visiblep.clear();

		H.erase(t2);
		H.insert(t);
		point tempr;
		facet tempt;
		for (int i = 0; i < 2; i++) {
			if (i) tempr = t.v1;
			else tempr = t.v2;
			if (r.pivot == tempr.pivot) {
				ProcessRidge(t1, r, t);
			}
			else {
				multimap<point, facet>::iterator it = M.find(tempr);
				if (it == M.end()) {
					M.insert(pair<point, facet>(tempr, t));
				}
				else {
					ProcessRidge(t, tempr, it->second);
				}
			}
		}
	}
	else {//右
		point p = (*t1iter->second.begin());
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
				ProcessRidge(t, r, t2);
			}
			else {
				multimap<point, facet>::iterator it = M.find(tempr);
				if (it == M.end()) {
					M.insert(pair<point, facet>(tempr, t));
				}
				else {
					ProcessRidge(it->second, tempr, t);
				}
			}
		}
	}
	return;
}