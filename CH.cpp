#include <iostream>
#include <cstdio>
#include <stdlib.h>
#include <vector>
#include <list>
#include <map>
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

bool *C = new bool[10000000000]();//


struct point
{
	double x, y;
}s[10000000];

typedef struct facet
{
	point v1, v2;
}Facet;


void RandmOnCircle(point *p, int n) {
	for (int i = 0; i < n; i++) {
		p[i].x = (double)n + n*cos((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
		p[i].y = (double)n + n*sin((((hash1(i)) % (n * 100))*360.0) / (100 * n)*PI / 180);
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
int cmp1(point a, point b) {
	if (a.y == b.y)
		return a.x < b.x;
		return a.y < b.y;
}

int cmp2(point a, point b) {
	if (atan2(a.y - s[0].y, a.x - s[0].x) != atan2(b.y - s[0].y, b.x - s[0].x))
			return (atan2(b.y - s[0].y, b.x - s[0].x)-(atan2(a.y - s[0].y, a.x - s[0].x)) > EPSINON);
		return (b.x - a.x > EPSINON) ;
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
		}
	}
	else if (type_of_input == 1) {
		RandmOnCircle(p, n);
	}
//	for (int i = 0; i < n; i++) cout << p[i].x << "," << p[i].y << endl;
	vector<Facet> H;
	while (collinear(p[0], p[1], p[2])) {
		int x = rand() % n;
		point temp = p[2];
		p[2] = p[x];
		p[x] = temp;
	}
	if (start(p[0], p[1], p[2])) {
		facet t = { p[0],p[2] };
		H.push_back(t);
		t = { p[2],p[1] };
		H.push_back(t);
		t = { p[1],p[0] };
		H.push_back(t);
	}
	else {
		facet t = { p[0],p[1] };
		H.push_back(t);
		t = { p[1],p[2] };
		H.push_back(t);
		t = { p[2],p[0] };
		H.push_back(t);
	}
	cout << "check" << endl;
	cout << sizeof(C[0]) << " " << sizeof(int) << endl;
	timer t1,t2,t3,t4,t5;
	
	for (int i = 0; i < 3; i++) {
		for (int j = 3; j < n; j++) {
			C[i*n + j - 3] = visible(p[j], H[i]);
		}
	}
	
	vector<Facet> R;
	int ccount = 3;
	int* tool = new int[n]();
	int location; 
	bool flag1, flag2, inside;
	facet lt, rt;
	int lr, rr;
	int temp;
	
	
	for (int i = 3; i < n; i++) {
		t1.start();
		//cout << i << endl;

		inside=false;
/*		for (int i = 0; i < n; i++){
         	       for(int j = 0; j<n;j++){
                        cout << C[i*n+j] << " ";
                	}
                cout << endl;
        	}		
*/		temp = H.size();
			if (C[(temp-1)*n + i - 3]) {//end==1
				for (int k = 0; k < temp; k++) {
					if (!C[k*n + i - 3]) {
						rr = k - 1;
						if(k == 0){
							rr = temp-1;
						}
						break;
					}
				}
				for (int k = 0; k < temp; k++) {
					if (!C[(temp - k - 1) * n + i - 3]) {
						lr = temp - k;
						break;
					}
				}
			}
			else {//end==0
				for (int k = 0; k < temp; k++) {
					if (C[k*n + i - 3]) {
						lr = k;
						break;
					}
					if (k == temp - 1)inside = true;
				}
				for (int k = 0; k < temp; k++) {
					if (C[(temp - 1 - k)*n + i - 3]) {
						rr = temp - k - 1;
						break;
					}
				}
			}
		t1.stop();
		if (inside) continue;
	t3.start();
	rt = { p[i], H[rr].v2 };
		lt = { H[lr].v1, p[i] };
	/*	cout << "rr:" << rr << endl;
		cout << "lr:" << lr << endl;
		cout << "rt" << rt.v1.x << "," << rt.v1.y<< " " << rt.v2.x << "," << rt.v2.y << endl;
		cout << "lt" << lt.v1.x << "," << lt.v1.y << " " << lt.v2.x << "," << lt.v2.y << endl;
*/
		if(lr>rr){
			H.erase(H.begin() + lr, H.end());
			H.erase(H.begin(), H.begin() + rr);
			H.push_back(lt);
                   	H[0] = rt;
			
			for (int j = (rr)*n; j < (lr+1) * n ; j++) {//向前移位
                                C[j-n*(rr)] = C[j];
                        }
			for (int j = i-3; j < n; j++) {
                                if (C[j] || C[n+j])
                                        C[j] = visible(p[j+3], rt);
                                if (C[(lr-rr-1)*n + j] || C[(lr-rr) * n + j])
                                        C[(lr-rr)*n+j] = visible(p[j+3], lt);
                        }
			ccount = lr-rr+1;

		}
		else{



			if (rr-lr == 0) {
				
				H[lr] = rt;
				H.insert(H.begin()+lr,lt);t5.start();
				//cout << ccount << endl;
				for (int j = ccount * n - 1; j > lr*n - 1; j--) {//向后移位
					C[j + n] = C[j];
				}t5.stop();
				for (int j = i-3; j < n; j++) {
					if (C[(lr+1)*n + j] || C[((lr+ccount)%(ccount+1))*n + j])
						C[lr*n + j] = visible(p[j+3], lt);
					if (C[(rr+1)*n + j] || C[((rr+2)%(ccount+1)) * n + j])
						//cout << "test" << p[j].x << "," << p[j].y << endl;
						C[(lr + 1)*n + j] = visible(p[j+3], rt);
				}
				ccount++; 
			}
			else if(rr - lr == 1){t4.start();
				H[lr] = lt;
				H[lr+1] = rt;
				for (int j = i-3; j < n; j++) {
					if (C[lr*n + j] || C[((lr-1+ccount)%(ccount)) * n + j])
						C[lr*n + j] = visible(p[j+3], lt);
					if (C[((rr+1)%(ccount))*n + j] || C[rr * n + j])
						C[(rr)*n + j] = visible(p[j+3], rt);
				}t4.stop();
			}
			
			else {
				H[lr] = lt;
				H[lr + 1] = rt;
				H.erase(H.begin() + lr + 2, H.begin() + rr+1);
				for (int j = i-3; j < n; j++) {
					if (C[lr*n + j] || C[((lr - 1+ccount)%(ccount)) * n + j])
						C[lr*n + j] = visible(p[j+3], lt);
					if (C[(rr)*n + j] || C[(rr + 1) * n + j])
						C[(rr)*n + j] = visible(p[j+3], rt);
				}
	//			cout << "check" << endl;
				for (int j = rr*n; j < ccount * n ; j++) {//向前移位
					C[j-n*(rr-lr-1)] = C[j];
				}
				ccount = ccount - rr + lr + 1;
			}
			
		} 
//for (int x = 0; x < H.size(); x++) cout << i << ": " << H[x].v1.x << "," << H[x].v1.y << " " << H[x].v2.x << "," << H[x].v2.y << endl;


	t3.stop();
	}
	
	t2.start();
			memset(s,0,sizeof(s));
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
		t2.stop();
		sort(s,s+top+1,cmp1);
//		for (int i = 0; i < top+1; i++) cout << s[i].x << "," << s[i].y << endl;
		for(int i = 0; i < top; i++){
//			cout << i << ":" << top << endl;
			while(s[i].x == s[i+1].x && s[i].y == s[i+1].y){
				for(int j = i+1; j < top; j++){
					s[j] = s[j+1];
				}
				top--;
			}			
		}
		if(H.size() >= top+1){	
			for(int i = 0; i < top+1; i++){
				for(int j = 0; j < H.size(); j++){
					if(s[i].x == H[j].v1.x && s[i].y == H[j].v1.y) break;
					if(j == H.size()-1)return 0;
				}
			}
			cout << "correct" << endl;
		}/**/
//	for (int i = 0; i < H.size(); i++) cout << H[i].v1.x << " " << H[i].v1.y << " " << H[i].v2.x << " " << H[i].v2.y << endl;
//for (int i = 0; i < top+1; i++) cout << s[i].x << "," << s[i].y << endl;


		if (H.size() >= top + 1) {
			cout << "size1: " << H.size() << endl;
			cout << "size2: " << top + 1 << endl;
			cout << "time1: " << t1.get_total() << endl;
			cout << "time2: " << t2.get_total() << endl;
			cout << "time3: " << t3.get_total() << endl;
			cout << "time4: " << t4.get_total() << endl;
			cout << "time5: " << t5.get_total() << endl;
		}
		else {
			for (int i = 0; i < top + 1; i++) {
				for (int j = 0; j < H.size(); j++) {
					if (s[i].x == H[j].v1.x && s[i].y == H[j].v1.y) break;
					if (j == H.size() - 1)cout << s[i].x << "," << s[i].y << endl;;
				}
			}
		}
	

	return 0;
}
