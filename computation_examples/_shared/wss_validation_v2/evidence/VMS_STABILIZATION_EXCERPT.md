# 实际冻结流体源码：三维VMS中的时间步与连续性项

来源：/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular/wss_audit/evidence/source/solver_frozen/fluid.cpp

SHA256 `df04669e8393de6f1eef8fcc18ce883dd39369bad64c3888f665b5525ef2ee64`。这是原冻结源码的只读摘录，不在本轮修改求解器。

```cpp
1490:   double ux[3][3] = {};
1491:   double uxx[3][3][3] = {};
1492: 
1493:   for (int a = 0; a < eNoNw; a++) {
1494:     ud[0] = ud[0] + Nw(a)*(al(0,a)-bfl(0,a));
1495:     ud[1] = ud[1] + Nw(a)*(al(1,a)-bfl(1,a));
1496:     ud[2] = ud[2] + Nw(a)*(al(2,a)-bfl(2,a));
1497: 
1498:     u[0] = u[0] + Nw(a)*yl(0,a);
1499:     u[1] = u[1] + Nw(a)*yl(1,a);
1500:     u[2] = u[2] + Nw(a)*yl(2,a);
1501: 
1502:     ux[0][0] = ux[0][0] + Nwx(0,a)*yl(0,a);
1503:     ux[1][0] = ux[1][0] + Nwx(1,a)*yl(0,a);
1504:     ux[2][0] = ux[2][0] + Nwx(2,a)*yl(0,a);
1505:     ux[0][1] = ux[0][1] + Nwx(0,a)*yl(1,a);
1506:     ux[1][1] = ux[1][1] + Nwx(1,a)*yl(1,a);
1507:     ux[2][1] = ux[2][1] + Nwx(2,a)*yl(1,a);
1508:     ux[0][2] = ux[0][2] + Nwx(0,a)*yl(2,a);
1509:     ux[1][2] = ux[1][2] + Nwx(1,a)*yl(2,a);
1510:     ux[2][2] = ux[2][2] + Nwx(2,a)*yl(2,a);
1511: 
1512:     uxx[0][0][0] += Nwxx(0,a)*yl(0,a);
1513:     uxx[1][0][1] += Nwxx(1,a)*yl(0,a);
1514:     uxx[2][0][2] += Nwxx(2,a)*yl(0,a);
1515:     uxx[1][0][0] += Nwxx(3,a)*yl(0,a);
```

```cpp
1648:   //
1649:   double up[3] = {};
1650:   double updu[3][3][MAX_SIZE] = {};
1651:   double tauM = 0.0;
1652: 
1653:   if (vmsFlag) {
1654:     // Stabilization parameters
1655:     double kT = 4.0 * pow(ctM/dt,2.0);
1656:     
1657:     // If we consider the NSB model, we need to add an extra term inside the computation for the stab parameter 
1658:     kT = kT + pow(brinkman_inverse_permeability*mu/rho, 2.0);
1659:     
1660:     // In case of unfitted RIS, compute the delta function at the quad point,
1661:     // add the additional value to the stabilization param 
1662:     kT = kT + pow(urisFactorTotal, 2.0);
1663: 
1664:     double kU = u[0]*u[0]*Kxi(0,0) + u[1]*u[0]*Kxi(1,0) + u[2]*u[0]*Kxi(2,0)
1665:               + u[0]*u[1]*Kxi(0,1) + u[1]*u[1]*Kxi(1,1) + u[2]*u[1]*Kxi(2,1)
1666:               + u[0]*u[2]*Kxi(0,2) + u[1]*u[2]*Kxi(1,2) + u[2]*u[2]*Kxi(2,2);
1667: 
1668:     double kS = Kxi(0,0)*Kxi(0,0) + Kxi(1,0)*Kxi(1,0) + Kxi(2,0)*Kxi(2,0)
1669:               + Kxi(0,1)*Kxi(0,1) + Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1)
1670:               + Kxi(0,2)*Kxi(0,2) + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2);
1671: 
1672:     kS = ctC * kS * pow(mu/rho,2.0);
1673:     tauM = 1.0 / (rho * sqrt( kT + kU + kS ));
1674: 
1675:     double rV[3];
1676:     rV[0] = ud[0] + u[0]*ux[0][0] + u[1]*ux[1][0] + u[2]*ux[2][0];
1677:     rV[1] = ud[1] + u[0]*ux[0][1] + u[1]*ux[1][1] + u[2]*ux[2][1];
1678:     rV[2] = ud[2] + u[0]*ux[0][2] + u[1]*ux[1][2] + u[2]*ux[2][2];
1679: 
1680:     double rS[3];
1681:     rS[0] = mu_x[0]*es[0][0] + mu_x[1]*es[1][0] + mu_x[2]*es[2][0] + mu*d2u2[0];
1682:     rS[1] = mu_x[0]*es[0][1] + mu_x[1]*es[1][1] + mu_x[2]*es[2][1] + mu*d2u2[1];
1683:     rS[2] = mu_x[0]*es[0][2] + mu_x[1]*es[1][2] + mu_x[2]*es[2][2] + mu*d2u2[2];
1684: 
1685:     // up[0] = -tauM*(rho*rV[0] + px[0] - rS[0] + mu*brinkman_inverse_permeability*u[0]);
1686:     // up[1] = -tauM*(rho*rV[1] + px[1] - rS[1] + mu*brinkman_inverse_permeability*u[1]);
1687:     // up[2] = -tauM*(rho*rV[2] + px[2] - rS[2] + mu*brinkman_inverse_permeability*u[2]);
1688: 
1689:     up[0] = -tauM*(rho*rV[0] + px[0] - rS[0] + mu*brinkman_inverse_permeability*u[0]
1690:                    + urisFactorTotal*u[0] - urisValveVelTermTotal[0]);
1691:     up[1] = -tauM*(rho*rV[1] + px[1] - rS[1] + mu*brinkman_inverse_permeability*u[1]
1692:                    + urisFactorTotal*u[1] - urisValveVelTermTotal[1]);
1693:     up[2] = -tauM*(rho*rV[2] + px[2] - rS[2] + mu*brinkman_inverse_permeability*u[2]
1694:                    + urisFactorTotal*u[2] - urisValveVelTermTotal[2]);
1695: 
1696:     for (int a = 0; a < eNoNw; a++) {
1697:       double uNx = u[0]*Nwx(0,a) + u[1]*Nwx(1,a) + u[2]*Nwx(2,a);
1698:       // T1 = -rho*uNx + mu*(Nwxx(0,a) + Nwxx(1,a) + Nwxx(2,a)) + mu_x[0]*Nwx(0,a) + mu_x[1]*Nwx(1,a) + mu_x[2]*Nwx(2,a) - mu*brinkman_inverse_permeability*Nw(a);
1699: 
1700:       T1 = -rho*uNx + mu*(Nwxx(0,a) + Nwxx(1,a) + Nwxx(2,a)) 
1701:            + mu_x[0]*Nwx(0,a) + mu_x[1]*Nwx(1,a) + mu_x[2]*Nwx(2,a) 
1702:            - mu*brinkman_inverse_permeability*Nw(a)
1703:            - urisFactorTotal*Nw(a);
1704: 
1705:       updu[0][0][a] = mu_x[0]*Nwx(0,a) + d2u2[0]*mu_g*esNx[0][a] + T1;
1706:       updu[1][0][a] = mu_x[1]*Nwx(0,a) + d2u2[1]*mu_g*esNx[0][a];
1707:       updu[2][0][a] = mu_x[2]*Nwx(0,a) + d2u2[2]*mu_g*esNx[0][a];
1708:   
1709:       updu[0][1][a] = mu_x[0]*Nwx(1,a) + d2u2[0]*mu_g*esNx[1][a];
1710:       updu[1][1][a] = mu_x[1]*Nwx(1,a) + d2u2[1]*mu_g*esNx[1][a] + T1;
1711:       updu[2][1][a] = mu_x[2]*Nwx(1,a) + d2u2[2]*mu_g*esNx[1][a];
1712:   
1713:       updu[0][2][a] = mu_x[0]*Nwx(2,a) + d2u2[0]*mu_g*esNx[2][a];
1714:       updu[1][2][a] = mu_x[1]*Nwx(2,a) + d2u2[1]*mu_g*esNx[2][a];
1715:       updu[2][2][a] = mu_x[2]*Nwx(2,a) + d2u2[2]*mu_g*esNx[2][a] + T1;
1716:     }
1717: 
1718:   } else {
1719:     tauM = 0.0;
1720:     std::memset(up, 0, sizeof up);
1721:     std::memset(updu, 0, sizeof updu);
1722:   }
1723: 
1724:   //  Local residual
1725:   //
```

```cpp
1995:   }
1996:   std::transform(mu_x.begin(), mu_x.end(), mu_x.begin(), [mu_g](double &v){return mu_g*v;});
1997:   //mu_x(:) = mu_g * mu_x(:)
1998: 
1999:   // Stabilization parameters
2000:   //
2001:   double kT = 4.0 * pow(ctM/dt,2.0);
2002:   
2003:   // If we consider the NSB model, we need to add an extra term inside the computation for the stab parameter 
2004:   kT = kT + pow(brinkman_inverse_permeability*mu/rho, 2.0);
2005: 
2006:   // In case of unfitted RIS, compute the delta function at the quad point,
2007:   // add the additional value to the stabilization param 
2008:   kT = kT + pow(urisFactorTotal, 2.0);
2009: 
2010:   double kU = u[0]*u[0]*Kxi(0,0) + u[1]*u[0]*Kxi(1,0) + u[2]*u[0]*Kxi(2,0)
2011:             + u[0]*u[1]*Kxi(0,1) + u[1]*u[1]*Kxi(1,1) + u[2]*u[1]*Kxi(2,1)
2012:             + u[0]*u[2]*Kxi(0,2) + u[1]*u[2]*Kxi(1,2) + u[2]*u[2]*Kxi(2,2);
2013: 
2014:   double kS = Kxi(0,0)*Kxi(0,0) + Kxi(1,0)*Kxi(1,0) + Kxi(2,0)*Kxi(2,0)
2015:             + Kxi(0,1)*Kxi(0,1) + Kxi(1,1)*Kxi(1,1) + Kxi(2,1)*Kxi(2,1)
2016:             + Kxi(0,2)*Kxi(0,2) + Kxi(1,2)*Kxi(1,2) + Kxi(2,2)*Kxi(2,2);
2017:   kS = ctC * kS * pow(mu/rho,2.0);
2018:   double tauM = 1.0 / (rho * sqrt( kT + kU + kS ));
2019: 
2020:   #ifdef debug_fluid_3d_m
2021:   dmsg << "kT: " << kT;
2022:   dmsg << "kU: " << kU;
2023:   dmsg << "kS: " << kS;
2024:   dmsg << "tauM: " << tauM;
2025:   #endif
2026: 
2027:   double rV[3] = {};
2028:   rV[0] = ud[0] + u[0]*ux[0][0] + u[1]*ux[1][0] + u[2]*ux[2][0];
2029:   rV[1] = ud[1] + u[0]*ux[0][1] + u[1]*ux[1][1] + u[2]*ux[2][1];
2030:   rV[2] = ud[2] + u[0]*ux[0][2] + u[1]*ux[1][2] + u[2]*ux[2][2];
2031: 
2032:   double rS[3] = {};
2033:   rS[0] = mu_x[0]*es[0][0] + mu_x[1]*es[1][0] + mu_x[2]*es[2][0] + mu*d2u2[0];
2034:   rS[1] = mu_x[0]*es[0][1] + mu_x[1]*es[1][1] + mu_x[2]*es[2][1] + mu*d2u2[1];
2035:   rS[2] = mu_x[0]*es[0][2] + mu_x[1]*es[1][2] + mu_x[2]*es[2][2] + mu*d2u2[2];
2036: 
2037:   double up[3] = {};
2038:   // up[0] = -tauM*(rho*rV[0] + px[0] - rS[0] + mu*brinkman_inverse_permeability * u[0]);
2039:   // up[1] = -tauM*(rho*rV[1] + px[1] - rS[1] + mu*brinkman_inverse_permeability * u[1]);
2040:   // up[2] = -tauM*(rho*rV[2] + px[2] - rS[2] + mu*brinkman_inverse_permeability * u[2]);
2041: 
2042:   up[0] = -tauM*(rho*rV[0] + px[0] - rS[0] + mu*brinkman_inverse_permeability * u[0]
2043:                  + urisFactorTotal * u[0] - urisValveVelTermTotal[0]);
2044:   up[1] = -tauM*(rho*rV[1] + px[1] - rS[1] + mu*brinkman_inverse_permeability * u[1]
2045:                  + urisFactorTotal * u[1] - urisValveVelTermTotal[1]);
2046:   up[2] = -tauM*(rho*rV[2] + px[2] - rS[2] + mu*brinkman_inverse_permeability * u[2]
2047:                  + urisFactorTotal * u[2] - urisValveVelTermTotal[2]);
2048: 
2049:   double tauC, tauB, pa;
2050:   double eps = std::numeric_limits<double>::epsilon();
2051:   double ua[3] = {};
2052: 
2053:   if (vmsFlag) {
2054:     tauC = 1.0 / (tauM * (Kxi(0,0) + Kxi(1,1) + Kxi(2,2)));
2055:     tauB = up[0]*up[0]*Kxi(0,0) + up[1]*up[0]*Kxi(1,0)
2056:          + up[2]*up[0]*Kxi(2,0) + up[0]*up[1]*Kxi(0,1)
2057:          + up[1]*up[1]*Kxi(1,1) + up[2]*up[1]*Kxi(2,1)
2058:          + up[0]*up[2]*Kxi(0,2) + up[1]*up[2]*Kxi(1,2)
2059:          + up[2]*up[2]*Kxi(2,2);
2060: 
2061:     if (utils::is_zero(tauB)) {
2062:       tauB = eps;
2063:     }
2064:     tauB = rho / sqrt(tauB);
2065: 
```

## 连续性弱残量（3D）

```cpp
1725:   //
1726:   for (int a = 0; a < eNoNq; a++) {
1727:     double upNx = up[0]*Nqx(0,a) + up[1]*Nqx(1,a) + up[2]*Nqx(2,a);
1728:     lR(3,a) = lR(3,a) + w*(Nq(a)*divU - upNx);
1729:   }
1730: 
```

此残量含稳定化细尺度通量项，收敛不等同于原P1速度逐点div(u)=0。控制体积通量/散度检查直接测量实际P1场，不用代数残差替代它。
## 实际调用选择与来源核对

服务器构建源 fluid.cpp、本地 external/flow_solver_source 及上一轮冻结副本SHA256均为`df04669e8393de6f1eef8fcc18ce883dd39369bad64c3888f665b5525ef2ee64`。当前XML没有启用Taylor-Hood；实际单元TET4，速度/压力等阶。

### fluid.cpp

SHA256 `df04669e8393de6f1eef8fcc18ce883dd39369bad64c3888f665b5525ef2ee64`

```cpp
491:   using namespace consts;
492: 
493:   const int eNoN = lM.eNoN;
494:   bool vmsStab = false;
495: 
496:   if (lM.nFs == 1) {
497:      vmsStab = true;
498:   } else {
499:      vmsStab = false;
500:   }
501: 
```

### fluid.cpp

SHA256 `df04669e8393de6f1eef8fcc18ce883dd39369bad64c3888f665b5525ef2ee64`

```cpp
725:       // Compute continuity residual and tangent matrix.
726:       //
727:       if (nsd == 3) {
728:         auto N0 = fs[0].N.rcol(g); 
729:         auto N1 = fs[1].N.rcol(g); 
730:         double urisFactorTotal = 0.0;
731:         Vector<double> urisValveVelTermTotal(nsd);
732:         if (com_mod.urisFlag) {
733:           urisFactorTotal = urisFactorTotalEl(g);
734:           urisValveVelTermTotal = urisValveVelTermTotalEl.rcol(g);
735:         }
736:         fluid_3d_c(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, 
737:               Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK, brinkman_inverse_permeability,
738:               urisFactorTotal, urisValveVelTermTotal);
```

### read_files.cpp

SHA256 `57ab378dc3f7e8442b28b99a908b91a546db2f63a769244a1a9e2510d7e624b4`

```cpp
1740:   bool THflag = false; 
1741:   if (eq_params->use_taylor_hood_type_basis.defined()) { 
1742:     THflag = eq_params->use_taylor_hood_type_basis.value(); 
1743:   }
1744:   EquationProps propL{consts::PhysicalPropertyType::NA};
1745:   EquationOutputs outPuts;
1746:   EquationNdop nDOP;
1747: 
1748:   set_equation_properties(simulation, eq_params, lEq, propL, outPuts, nDOP);
1749: 
1750:   // Check if explicit geometric coupling is allowed for the equation
1751:   if (lEq.expl_geom_cpl) {
1752:     if (lEq.phys != EquationType::phys_FSI) {
1753:       throw std::runtime_error("Explicit geometric coupling is only allowed for FSI equation.");
1754:     }
1755:   }
1756: 
1757:   // Read parameters related to VTU output.
1758:   read_outputs(simulation, eq_params, lEq, nDOP, outPuts);
1759: 
1760:   // Set the number of function spaces
1761:   for (int iM = 0; iM < com_mod.nMsh; iM++) {
1762:     auto& msh = com_mod.msh[iM];
1763:     if (!THflag) {
1764:       msh.nFs = 1;
1765:     } else { 
1766:       if (com_mod.ibFlag) {
1767:         throw std::runtime_error("Taylor-Hood basis is not implemented for immersed boundaries.");
1768:       }
1769:       if (msh.lShl || msh.lFib || (msh.eType == ElementType::NRB)) {
```

