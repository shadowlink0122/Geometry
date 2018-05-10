
#include <cmath>
#include <iomanip>
using namespace std;

#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

namespace Geometry{
  #define equals(a,b) (fabs((a) - (b)) < EPS)
  #define BOTTOM 0
  #define LEFT 1
  #define RIGHT 2
  #define TOP 3
  #define gprint(X) cout << fixed << setprecision(0) << (X) << endl
  #define gprints(X,Y) cout << fixed << setprecision(0) << (X) << " " << (Y) << endl
  const double EPS = 1e-10;
  static const int COUNTER_CLOCKWISE = 1;
  static const int CLOCKWISE = -1;
  static const int ONLINE_BACK = 2;
  static const int ONLINE_FRONT = -2;
  static const int ON_SEGMENT = 0;
  
  class Point{
  public:
    double x,y;
    Point(double x = 0,double y = 0): x(x),y(y){}
    Point operator + (Point p){ return Point(x + p.x , y + p.y); }
    Point operator - (Point p){ return Point(x - p.x , y - p.y); }
    Point operator += (Point p){ return Point(x += p.x,y += p.y); }
    Point operator -= (Point p){ return Point(x -= p.x,y -= p.y); }
    Point operator * (double a){ return Point(a * x , a * y); }
    Point operator / (double a){ return Point(x / a , y / a); }
    bool operator < (const Point &p) const{ return x != p.x ? x < p.x : y < p.y; }
    bool operator == (const Point &p) const{ return equals(x,p.x) && equals(y,p.y); }

    double abs(){ return sqrt(norm()); }
    double norm(){ return x * x + y * y; }
    friend ostream& operator << (ostream &os,Point &p){ return os << p.x << " " << p.y; }
    friend istream& operator >> (istream &is,Point &p){ return is >> p.x >> p.y; }
  };

  struct Segment{
    Point p1,p2;
    Segment(){};
    Segment(Point p1,Point p2): p1(p1) , p2(p2){};
    friend ostream& operator << (ostream &os,Segment &seg){ return os << seg.p1 << " " << seg.p2; }
    friend istream& operator >> (istream &is,Segment &seg){ return is >> seg.p1 >> seg.p2; }
  };

  class Circle{
  public:
    Point c;
    double r;
    Circle(Point c = Point(),double r = 0.0): c(c),r(r){}
    friend ostream& operator << (ostream &os,Circle &c){ return os << c.c << " " << c.r; }
    friend istream& operator >> (istream &is,Circle &c){ return is >> c.c >> c.r; }
  };

  typedef Point Vector;
  typedef Segment Line;
  typedef vector<Point> Polygon;

  double norm(Vector a){ return  a.x * a.x + a.y * a.y; }
  double abs(Vector a){ return sqrt(norm(a)); }
  double dot(Vector a,Vector b){ return a.x * b.x + a.y * b.y; }
  double cross(Vector a,Vector b){ return a.x * b.y - a.y * b.x; }
  bool Orthogonal(Vector a,Vector b){ return equals(dot(a,b),0.0); }
  bool Orthogonal(Point a1,Point a2,Point b1,Point b2){ return Orthogonal(a1 - a2,b1 - b2); }
  bool Orthogonal(Segment s1,Segment s2){ return equals(dot(s1.p2 - s1.p1,s2.p2 - s2.p1) , 0.0); }
  bool Parallel(Vector a,Vector b){ return equals(cross(a,b),0.0); }
  bool Parallel(Point a1,Point a2,Point b1,Point b2){ return Parallel(a1 - a2 , b1 - b2); }
  bool Parallel(Segment s1,Segment s2){ return equals(cross(s1.p2 - s1.p1 , s2.p2 - s2.p1),0.0); }
  bool less(const Point &a,const Point &b){ return (a.x != b.x) ? a.x < b.x : a.y < b.y; }
  bool greater(const Point &a,const Point b){ return !less(a,b); }

  Point project(Segment s,Point p){
    Vector base = s.p2 - s.p1;
    double r = dot(p - s.p1 , base) / norm(base);
    return s.p1 + base * r;
  }

  Point reflect(Segment s,Point p){
    return p + (project(s,p) - p) * 2.0;
  }

  int ccw(Point &p0,Point p1,Point p2){
    p1 -= p0;
    p2 -= p0;
    if(cross(p1,p2) > EPS)return COUNTER_CLOCKWISE;
    if(cross(p1,p2) < -EPS)return CLOCKWISE;
    if(dot(p1,p2) < -EPS)return ONLINE_BACK;
    if(p1.norm() < p2.norm())return ONLINE_FRONT;
    return ON_SEGMENT;
  }

  bool intersect(Point p1,Point p2,Point p3,Point p4){
    return (ccw(p1,p2,p3) * ccw(p1,p2,p4) <= 0 && 
            ccw(p3,p4,p1) * ccw(p3,p4,p2) <= 0);
  }
  bool intersect(Segment s1,Segment s2){ return intersect(s1.p1 , s1.p2 , s2.p1 , s2.p2); }

  double getDistance(Vector a,Vector b){ return abs(a - b); }
  double getDistanceLP(Line l,Point p){ return abs(cross(l.p2 - l.p1 , p - l.p1) / abs(l.p2 - l.p1)); }
  double getDistanceSP(Segment s,Point p){
    if(dot(s.p2 - s.p1,p - s.p1) < 0.0) return abs(p - s.p1);
    if(dot(s.p1 - s.p2,p - s.p2) < 0.0) return abs(p - s.p2);
    return getDistanceLP(s,p);
  }
  double getDistance(Segment s1,Segment s2){
    if(intersect(s1,s2))return 0.0;
    return min(min(getDistanceSP(s1,s2.p1),getDistanceSP(s1,s2.p2)),
               min(getDistanceSP(s2,s1.p1),getDistanceSP(s2,s1.p2)));
  }

  Point getCrossPoint(Segment s1,Segment s2){
    Vector base = s2.p2 - s2.p1;
    double d1 = abs(cross(base,s1.p1 - s2.p1));
    double d2 = abs(cross(base,s1.p2 - s2.p1));
    double t = d1 / (d1 + d2);
    return s1.p1 + (s1.p2 - s1.p1) * t;
  }

  pair<Point,Point> getCrossPoints(Circle c,Line l){
    Vector pr = project(l,c.c);
    Vector e = (l.p2 - l.p1) / abs(l.p2 - l.p1);
    double base = sqrt(c.r * c.r - norm(pr - c.c));
    return make_pair(pr + e * base , pr - e * base);
  }

  double arg(Vector p){ return atan2(p.y,p.x); }
  Vector polar(double a,double r){ return Point(cos(r) * a, sin(r) * a); }

  pair<Point,Point> getCrossPoints(Circle c1,Circle c2){
    double d = abs(c1.c - c2.c);
    double a = acos((c1.r * c1.r + d * d - c2.r * c2.r) / (2 * c1.r * d));
    double t = arg(c2.c - c1.c);
    return make_pair(c1.c + polar(c1.r,t + a),c1.c + polar(c1.r,t - a));
  }
  
  int contains(Polygon g,Point p){
    int n = g.size();
    bool x = false;
    for(int i = 0;i < n;i++){
      Point a = g[i] - p, b = g[(i + 1) % n] - p;
      if(abs(cross(a,b)) < EPS && dot(a,b) < EPS)return 1;
      if(a.y > b.y)swap(a,b);
      if(a.y < EPS && EPS < b.y && cross(a,b) > EPS)x = !x;
    }
    return (x ? 2 : 0);
  }

  //=====================================================//
  //==================== Andrew Scan ====================//
  //=====================================================//
 
  int Accw(Point p0,Point p1,Point p2){
    Vector a = p1 - p0;
    Vector b = p2 - p0;
    if(cross(a,b) > EPS)return 1;
    return -1;
  }
 
  Polygon andrewScan(Polygon s){
    Polygon ans;
    Polygon u,l;
    if(s.size() < 3)return s;
    sort(s.begin(),s.end());
    u.push_back(s[0]);
    u.push_back(s[1]);
    l.push_back(s[s.size()-1]);
    l.push_back(s[s.size()-2]);
    for(int i = 2;i < s.size();i++){
      for(int n = u.size();n >= 2 && Accw(u[n-2],u[n-1],s[i]) != -1;n--){
        u.pop_back();
      }
      u.push_back(s[i]);
    }
 
    for(int i = s.size() - 3;i >= 0;i--){
      for(int n = l.size();n >= 2 && Accw(l[n-2],l[n-1],s[i]) != -1;n--){
        l.pop_back();
      }
      l.push_back(s[i]);
    }
    reverse(l.begin(),l.end());
    for(int i = u.size() - 2;i >= 1;i--){
      l.push_back(u[i]);
    }
    int f = 0;
    for(int i = 1;i < l.size();i++){
      if(l[f].y > l[i].y || (l[f].y == l[i].y && l[f].x > l[i].x)){
        f = i;
      }
    }
    for(int i = 0;i < l.size();i++){
      ans.push_back(Point((int)l[(i + f)%l.size()].x,(int)l[(i + f)%l.size()].y));
    }
    return ans;
  }

  //=====================================================//
  //======================== END ========================//
  //=====================================================//

  class EndPoint{
  public:
    Point p;
    int seg,st;
    EndPoint(){}
    EndPoint(Point p,int seg,int st): p(p),seg(seg),st(st){}

    bool operator < (const EndPoint &ep) const{
      if(p.y == ep.p.y){
        return st < ep.st;
      }else{
        return p.y < ep.p.y;
      }
    }
  };

  EndPoint EP[2 * 100000];

  int manhattanIntersection(vector<Segment> S){
    int n = S.size();
    for(int i = 0,k = 0;i < n;i++){
      if(S[i].p1.y == S[i].p2.y){
        if(S[i].p1.x > S[i].p2.x)swap(S[i].p1,S[i].p2);
      }else if(S[i].p1.y > S[i].p2.y){
        swap(S[i].p1,S[i].p2);
      }

      if(S[i].p1.y == S[i].p2.y){
        EP[k++] = EndPoint(S[i].p1,i,LEFT);
        EP[k++] = EndPoint(S[i].p2,i,RIGHT);
      }else{
        EP[k++] = EndPoint(S[i].p1,i,BOTTOM);
        EP[k++] = EndPoint(S[i].p2,i,TOP);
      }
    }

    sort(EP,EP + (2 * n));

    set<int> BT;
    BT.insert(1000000001);
    int cnt = 0;
    for(int i = 0;i < 2 * n;i++){
      if(EP[i].st == TOP){
        BT.erase(EP[i].p.x);
      }else if(EP[i].st == BOTTOM){
        BT.insert(EP[i].p.x);
      }else if(EP[i].st == LEFT){
        set<int>::iterator b = lower_bound(BT.begin(),BT.end(),S[EP[i].seg].p1.x);
        set<int>::iterator e = upper_bound(BT.begin(),BT.end(),S[EP[i].seg].p2.x);
        cnt += distance(b,e);
      }
    }
    return cnt;
  }
};

#endif //_GEOMETRY_H_
