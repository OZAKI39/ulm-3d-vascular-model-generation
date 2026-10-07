"""Float64 Warp kernels. No host sync, allocation, QP or bubble loop in Python."""
import warp as wp


@wp.func
def cmul(a: wp.vec2d, b: wp.vec2d):
    return wp.vec2d(a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0])


@wp.kernel
def field_response(x: wp.array(dtype=wp.vec3d), grid: wp.array2d(dtype=wp.vec2d),
                   origin: wp.vec3d, spacing: wp.vec3d, nx: int, ny: int, nz: int,
                   uniform: int, uniform_p: wp.array(dtype=wp.vec2d),
                   knots: wp.array(dtype=wp.float64), values: wp.array3d(dtype=wp.vec2d),
                   valid: wp.array(dtype=int), kv: wp.array2d(dtype=wp.vec2d), kval: wp.array(dtype=int),
                   harmonics: wp.array2d(dtype=wp.vec2d), force: wp.array(dtype=wp.vec3d),
                   pressure: wp.array(dtype=wp.vec2d), diag: wp.array2d(dtype=wp.float64),
                   error: wp.array(dtype=int)):
    i=wp.tid()
    p=wp.vec2d(wp.float64(0.0)); gx=p; gy=p; gz=p
    if uniform == 1:
        p=uniform_p[0]
    else:
        q=wp.cw_div(x[i]-origin,spacing)
        if q[0]<wp.float64(0.0) or q[1]<wp.float64(0.0) or q[2]<wp.float64(0.0) or q[0]>wp.float64(nx-1) or q[1]>wp.float64(ny-1) or q[2]>wp.float64(nz-1):
            wp.atomic_or(error,0,1)
        else:
            ix=wp.min(int(wp.floor(q[0])),nx-2); iy=wp.min(int(wp.floor(q[1])),ny-2); iz=wp.min(int(wp.floor(q[2])),nz-2)
            w=q-wp.vec3d(wp.float64(ix),wp.float64(iy),wp.float64(iz))
            for a in range(2):
                for b in range(2):
                    for c in range(2):
                        wa=wp.float64(1.0)-w[0]; wb=wp.float64(1.0)-w[1]; wc=wp.float64(1.0)-w[2]
                        if a==1: wa=w[0]
                        if b==1: wb=w[1]
                        if c==1: wc=w[2]
                        weight=wa*wb*wc
                        idx=((ix+a)*ny+iy+b)*nz+iz+c
                        p=p+weight*grid[idx,0]; gx=gx+weight*grid[idx,1]; gy=gy+weight*grid[idx,2]; gz=gz+weight*grid[idx,3]
    pressure[i]=p
    amp=wp.length(p)
    diag[i,0]=wp.min(diag[i,0],amp); diag[i,1]=wp.max(diag[i,1],amp)
    cell=int(-1); knot=int(-1)
    for k in range(knots.shape[0]):
        if amp==knots[k]: knot=k
        if k<knots.shape[0]-1:
            if amp>knots[k] and amp<knots[k+1]: cell=k
    ok=int(1)
    if amp>wp.float64(0.0):
        if knot>=0:
            ok=kval[knot]
        elif cell>=0:
            ok=valid[cell]
        else:
            ok=0
    if ok==0: wp.atomic_or(error,0,2)
    theta=wp.atan2(p[1],p[0])
    for h in range(8):
        v=wp.vec2d(wp.float64(0.0))
        if amp>wp.float64(0.0) and ok==1:
            if knot>=0:
                v=kv[knot,h]
            else:
                fraction=(amp-knots[cell])/(knots[cell+1]-knots[cell])
                v=(wp.float64(1.0)-fraction)*values[cell,0,h]+fraction*values[cell,1,h]
            angle=wp.float64(h+1)*theta
            v=cmul(v,wp.vec2d(wp.cos(angle),wp.sin(angle)))
        harmonics[i,h]=v
    v=harmonics[i,0]
    force[i]=-wp.float64(0.5)*wp.vec3d(wp.dot(v,gx),wp.dot(v,gy),wp.dot(v,gz))


@wp.kernel
def pair_force(x: wp.array(dtype=wp.vec3d), v: wp.array2d(dtype=wp.vec2d),
               force: wp.array(dtype=wp.vec3d), n: int, radius: wp.float64,
               omega: wp.float64, rho: wp.float64, regularize: int,
               transition: wp.float64, soft_core: wp.float64, error: wp.array(dtype=int)):
    b,i,j=wp.tid()
    if j>i:
        ii=b*n+i; jj=b*n+j
        delta=x[jj]-x[ii]; d2=wp.dot(delta,delta)
        if d2<=wp.float64(0.0):
            wp.atomic_or(error,0,4)
        else:
            if regularize==1:
                q=wp.min(wp.sqrt(d2)/(transition*wp.float64(2.0)*radius),wp.float64(1.0))
                core=soft_core*wp.float64(2.0)*radius*(wp.float64(1.0)-q)*(wp.float64(1.0)-q)
                d2=d2+core*core
            corr=wp.float64(0.0)
            for h in range(8):
                corr=corr+wp.float64((h+1)*(h+1))*wp.dot(v[ii,h],v[jj,h])
            c=rho*omega*omega*wp.float64(0.5)/(wp.float64(12.566370614359172)*d2*wp.sqrt(d2))*corr
            f=c*delta
            wp.atomic_add(force,ii,f); wp.atomic_add(force,jj,-f)


@wp.kernel
def free_velocity(force: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.vec3d), drag: wp.float64):
    i=wp.tid(); out[i]=force[i]/drag


@wp.func
def lub_phi(e: wp.float64):
    return wp.float64(0.5)/e-wp.float64(0.45)*wp.log(e)-wp.float64(3.0/56.0)*e*wp.log(e)


@wp.kernel
def local_mobility(x: wp.array(dtype=wp.vec3d), free: wp.array(dtype=wp.vec3d),
                   old: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.vec3d),
                   n: int, radius: wp.float64, activation: wp.float64, floor: wp.float64):
    i=wp.tid(); start=(i//n)*n
    mat=wp.identity(n=3,dtype=wp.float64); rhs=free[i]
    for j in range(start,start+n):
        if j!=i:
            delta=x[j]-x[i]; d=wp.length(delta)
            e=(d-wp.float64(2.0)*radius)/radius
            if e<activation and d>wp.float64(0.0):
                normal=delta/d
                c=wp.float64(0.5)*(lub_phi(wp.max(e,floor))-lub_phi(activation))
                mat=mat+c*wp.outer(normal,normal)
                rhs=rhs+c*wp.dot(normal,old[j])*normal
    out[i]=wp.float64(0.5)*(old[i]+wp.inverse(mat)*rhs)


@wp.kernel
def predict(x: wp.array(dtype=wp.vec3d), velocity: wp.array(dtype=wp.vec3d),
            out: wp.array(dtype=wp.vec3d), dt: wp.float64, diag: wp.array2d(dtype=wp.float64)):
    i=wp.tid(); dx=dt*velocity[i]; out[i]=x[i]+dx
    diag[i,5]=diag[i,5]+wp.length(dx)
    diag[i,7]=wp.max(diag[i,7],wp.length(velocity[i]))
    diag[i,9]=wp.max(diag[i,9],wp.length(dx))


@wp.kernel
def contact(x: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.vec3d), n: int,
            radius: wp.float64, relaxation: wp.float64, diag: wp.array2d(dtype=wp.float64),
            error: wp.array(dtype=int)):
    i=wp.tid(); start=(i//n)*n; dx=wp.vec3d(wp.float64(0.0)); penetration=wp.float64(0.0)
    for j in range(start,start+n):
        if i!=j:
            delta=x[i]-x[j]; d=wp.length(delta)
            overlap=wp.float64(2.0)*radius-d
            if overlap>wp.float64(0.0):
                penetration=wp.max(penetration,overlap)
                if d>wp.float64(0.0):
                    dx=dx+wp.float64(0.5)*overlap/d*delta
                else:
                    wp.atomic_or(error,0,4)
    dx=relaxation*dx
    out[i]=x[i]+dx
    diag[i,2]=wp.max(diag[i,2],penetration)
    diag[i,4]=diag[i,4]+wp.length(dx)
    diag[i,6]=wp.max(diag[i,6],wp.length(dx))


@wp.kernel
def accept(x: wp.array(dtype=wp.vec3d), out: wp.array(dtype=wp.vec3d), n: int,
           radius: wp.float64, diag: wp.array2d(dtype=wp.float64), error: wp.array(dtype=int)):
    i=wp.tid(); start=(i//n)*n; penetration=wp.float64(0.0)
    for j in range(start,start+n):
        if i!=j:
            penetration=wp.max(penetration,wp.float64(2.0)*radius-wp.length(x[i]-x[j]))
    if not wp.isfinite(x[i][0]) or not wp.isfinite(x[i][1]) or not wp.isfinite(x[i][2]):
        wp.atomic_or(error,0,8)
    diag[i,3]=wp.max(diag[i,3],penetration)
    out[i]=x[i]
