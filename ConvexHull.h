#pragma once
#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <atomic>
#include <mutex>
#include <string.h>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
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

unordered_map<facet, vector<point>> mapC;
unordered_map<point, facet> M;
unordered_set<Facet> H;

void init(point *p, int n){
	while (collinear(p[0], p[1], p[2])) {//排除初始三点共线
		int x = rand() % n;
		point temp = p[2];
		p[2] = p[x];
		p[x] = temp;
	}
	vector<point> visiblep;
	point temp;
	facet t;
	//确定每条线中的点顺时针排列
	//初始C,H
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
	
	//初始M
	/*facet tempf;
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
	}*/
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
		//cout << "3check3" << endl;
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
