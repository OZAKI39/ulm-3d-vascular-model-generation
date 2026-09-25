from pathlib import Path
S=Path(__file__).resolve().parents[1]
p=S/'src/wall_distance.hpp';s=p.read_text().replace('class WallDistance {','struct ClosestWallPoint { Vec point{}; double distance=0; bool ambiguous=false; };\nclass WallDistance {').replace(' double distance(Vec p) const;',' double distance(Vec p) const;\n ClosestWallPoint closest(Vec p) const;\n bool segment_safe_exact(Vec a,Vec b,double minimum) const;').replace(' int build(int begin,int end);',' void closest_search(int index,Vec p,ClosestWallPoint& best) const;\n bool segment_search(int index,Vec a,Vec b,double minimum) const;\n int build(int begin,int end);');p.write_text(s)
p=S/'src/wall_distance.cpp';s=p.read_text();s+=r'''
namespace passive {
static Vec closest_edge(Vec p,Vec a,Vec b){Vec u=minus(b,a);double den=dot(u,u);return plus_scaled(a,u,den>0?std::clamp(dot(minus(p,a),u)/den,0.,1.):0.);}
static Vec closest_triangle(Vec p,Triangle const&t){
 Vec ab=minus(t.b,t.a),ac=minus(t.c,t.a),ap=minus(p,t.a);double aa=dot(ab,ab),bb=dot(ab,ac),cc=dot(ac,ac),pa=dot(ap,ab),pc=dot(ap,ac),det=aa*cc-bb*bb;
 if(det>0){double v=(cc*pa-bb*pc)/det,w=(aa*pc-bb*pa)/det;if(v>=0&&w>=0&&v+w<=1)return plus_scaled(plus_scaled(t.a,ab,v),ac,w);}
 Vec best=closest_edge(p,t.a,t.b);for(auto q:{closest_edge(p,t.b,t.c),closest_edge(p,t.c,t.a)})if(norm(minus(p,q))<norm(minus(p,best)))best=q;return best;
}
void WallDistance::closest_search(int i,Vec p,ClosestWallPoint&best) const {
 auto const&n=nodes[i];double tol=64*std::numeric_limits<double>::epsilon()*std::max(norm(p),1e-6);
 if(box2(n,p)>(best.distance+tol)*(best.distance+tol))return;
 if(n.left<0){for(int k=n.begin;k<n.end;k++){Vec q=closest_triangle(p,triangles[k]);double d=norm(minus(p,q));
  if(d<best.distance-tol){best={q,d,false};}
  else if(std::abs(d-best.distance)<=tol){if(norm(minus(q,best.point))>8*tol)best.ambiguous=true;if(d<best.distance){best.point=q;best.distance=d;}}
 }return;}
 int a=n.left,b=n.right;if(box2(nodes[b],p)<box2(nodes[a],p))std::swap(a,b);closest_search(a,p,best);closest_search(b,p,best);
}
ClosestWallPoint WallDistance::closest(Vec p)const{ClosestWallPoint best;best.distance=INFINITY;closest_search(0,p,best);return best;}
static double segment_segment2(Vec p,Vec q,Vec a,Vec b){
 Vec u=minus(q,p),v=minus(b,a),w=minus(p,a);double A=dot(u,u),B=dot(u,v),C=dot(v,v),D=dot(u,w),E=dot(v,w),det=A*C-B*B;
 double best=std::min({edge2(p,a,b),edge2(q,a,b),edge2(a,p,q),edge2(b,p,q)});
 if(det>0){double s=(B*E-C*D)/det,t=(A*E-B*D)/det;if(s>=0&&s<=1&&t>=0&&t<=1){Vec d=minus(plus_scaled(p,u,s),plus_scaled(a,v,t));best=std::min(best,dot(d,d));}}
 return best;
}
static Vec cross(Vec a,Vec b){return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
static double segment_triangle2(Vec a,Vec b,Triangle const&t){
 Vec d=minus(b,a),n=cross(minus(t.b,t.a),minus(t.c,t.a));double den=dot(d,n);
 if(den!=0){double f=dot(minus(t.a,a),n)/den;if(f>=0&&f<=1){Vec x=plus_scaled(a,d,f),v=minus(t.b,t.a),w=minus(t.c,t.a),h=minus(x,t.a);double vv=dot(v,v),vw=dot(v,w),ww=dot(w,w),hv=dot(h,v),hw=dot(h,w),det=vv*ww-vw*vw;if(det>0){double u=(ww*hv-vw*hw)/det,z=(vv*hw-vw*hv)/det;if(u>=0&&z>=0&&u+z<=1)return 0;}}}
 return std::min({triangle2(a,t),triangle2(b,t),segment_segment2(a,b,t.a,t.b),segment_segment2(a,b,t.b,t.c),segment_segment2(a,b,t.c,t.a)});
}
bool WallDistance::segment_search(int i,Vec a,Vec b,double minimum)const{
 auto const&n=nodes[i];double lower2=0;for(int j=0;j<3;j++){double lo=std::min(a[j],b[j]),hi=std::max(a[j],b[j]),gap=std::max({n.lo[j]-hi,lo-n.hi[j],0.});lower2+=gap*gap;}if(lower2>=minimum*minimum)return true;
 if(n.left<0){for(int k=n.begin;k<n.end;k++)if(segment_triangle2(a,b,triangles[k])<minimum*minimum)return false;return true;}
 return segment_search(n.left,a,b,minimum)&&segment_search(n.right,a,b,minimum);
}
bool WallDistance::segment_safe_exact(Vec a,Vec b,double minimum)const{return segment_search(0,a,b,minimum);}
}
''';p.write_text(s)
