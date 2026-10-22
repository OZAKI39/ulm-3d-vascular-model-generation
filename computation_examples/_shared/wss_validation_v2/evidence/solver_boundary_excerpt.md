## external/flow_solver_source/Code/Source/solver/VtkData.cpp

SHA256: 888b1826e32d1a63318e0ed6f3f829e6c7586ab2b457766e79d5c9808eab711a

```cpp
240: void VtkData::set_points(const Array<double> &points) {
241:   const int num_coords = points.ncols();
242:   svmp::check<svmp::FE::InvalidArgumentException>(
243:       num_coords != 0, "The number of points is zero.");
244: 
245:   svmp::check<svmp::FE::InvalidArgumentException>(
246:       points.nrows() >= 3,
247:       "The point coordinates are given as an array of " +
248:           std::to_string(points.nrows()) +
249:           " rows, while three coordinates per point are needed.");
250: 
251:   auto node_coords = vtkSmartPointer<vtkPoints>::New();
252:   node_coords->Allocate(num_coords, 1000);
253:   node_coords->SetNumberOfPoints(num_coords);
254: 
255:   for (int i = 0; i < num_coords; i++) {
256:     node_coords->SetPoint(i, points(0, i), points(1, i), points(2, i));
257:   }
258: 
259:   vtk_data->SetPoints(node_coords);
260: }
261: 
262: void VtkData::set_connectivity(const int nsd, const Array<int> &conn) {
263:   int num_elems = conn.ncols();
264:   int np_elem = conn.nrows();
```

```cpp
330: void VtkData::copy_points(Array<double> &points) const {
331:   auto vtk_points = vtk_data->GetPoints();
332:   auto num_points = vtk_points->GetNumberOfPoints();
333: 
334:   svmp::check<svmp::FE::InvalidArgumentException>(
335:       points.nrows() >= 3 && points.ncols() >= num_points,
336:       "The " + std::to_string(num_points) + " points of the VTK file '" +
337:           file_name_ + "' do not fit in an array of " +
338:           std::to_string(points.nrows()) + " rows and " +
339:           std::to_string(points.ncols()) + " columns.");
340: 
341:   double point[3];
342:   for (int i = 0; i < num_points; i++) {
343:     vtk_points->GetPoint(i, point);
344:     points(0, i) = point[0];
345:     points(1, i) = point[1];
346:     points(2, i) = point[2];
347:   }
348: }
349: 
350: void VtkData::copy_point_data(const std::string &data_name,
351:                               Array<double> &mesh_data) const {
352:   const auto vtk_array = vtk_data->GetPointData()->GetArray(data_name.c_str());
353:   svmp::check<svmp::FE::InvalidArgumentException>(
```

## external/flow_solver_source/Code/Source/solver/read_files.cpp

SHA256: 57ab378dc3f7e8442b28b99a908b91a546db2f63a769244a1a9e2510d7e624b4

```cpp
140:   } else if (std::set<std::string>{"Neumann", "Neu"}.count(bc_type)) {
141:     lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_Neu)); 
142:     coupled_bc_type = BoundaryConditionType::bType_Neu;
143:     if ((lEq.phys == EquationType::phys_fluid) || (lEq.phys == EquationType::phys_FSI)) {
144:       lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_bfs));
145:     }
146: 
147:   } else if (std::set<std::string>{"Traction", "Trac"}.count(bc_type)) {
148:     lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_trac)); 
149: 
150:     if (bc_params->traction_values_file_path.defined()) { 
151:       int iM = lBc.iM;
152:       int iFa = lBc.iFa;
153:       lBc.gm.nTP = 2;
154:       auto& com_mod = simulation->com_mod;
155:       int nsd = com_mod.nsd;
156:       lBc.gm.dof = nsd;
157: 
158:       lBc.gm.t.resize(2); 
159:       lBc.gm.d.resize(nsd, com_mod.msh[iM].fa[iFa].nNo, 2);
160: 
161:       lBc.gm.t[0] = 0.0;
162:       lBc.gm.t[1] = 1.E+10;
163:       lBc.gm.period = lBc.gm.t[1];
164:       auto file_name = bc_params->traction_values_file_path.value();
165: 
166:       read_trac_bcff(com_mod, lBc.gm, com_mod.msh[iM].fa[iFa], file_name);
167: 
168:       lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_gen));
169:       lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_flat));
170: 
171:       if (bc_params->traction_multiplier.defined()) { 
172:         double rtmp = bc_params->traction_multiplier.value(); 
173:         lBc.gm.d *= rtmp;
174:       }
175: 
176:       lBc.eDrn.resize(nsd); 
177:       lBc.h.resize(nsd);
178:       lBc.weakDir = false;
179:       return; 
180:     }
181: 
182:   } else if (std::set<std::string>{"Robin", "Rbn"}.count(bc_type)) {
183:     lBc.bType = utils::ibset(lBc.bType, enum_int(BoundaryConditionType::bType_Robin)); 
```

```cpp
2459: void read_spatial_values(const ComMod& com_mod, const mshType& msh, const faceType& lFa, 
2460:     const std::string& file_name, bcType& lBc)
2461: {
2462:   std::ifstream file_stream;
2463:   file_stream.open(file_name);
2464:   if (!file_stream.is_open()) {
2465:     throw std::runtime_error("Failed to open the spatial values file '" + file_name + "'.");
2466:   }
2467: 
2468:   lBc.gx.resize(lFa.nNo); 
2469: 
2470:   // Preparing the pointer array
2471:   //
2472:   Vector<int> ptr(msh.gnNo);
2473:   ptr = -1;
2474:   for (int a = 0; a < lFa.nNo; a++) {
2475:     lBc.gx[a] = 0.0;
2476:     int Ac = lFa.gN[a];
2477:     Ac = msh.lN[Ac];
2478:     if (Ac == -1) {
2479:       throw std::runtime_error("Incorrect global node number detected for BC for mesh '" + msh.name + 
2480:           "' and face '" + lFa.name + "' for node " + std::to_string(a) + ".");
2481:     }
2482:     ptr[Ac] = a;
2483:   }
2484: 
2485:   for (int b = 0; b < lFa.nNo; b++) {
2486:     int Ac;
2487:     double rtmp;
2488:     file_stream >> Ac >> rtmp;
2489:     if ((Ac >= msh.gnNo) || (Ac < 0)) {
2490:       throw std::runtime_error("The node number " + std::to_string(Ac) + 
2491:             " in the spatial values file '" + file_name + " is larger than the number of nodes in the mesh.");
2492:     }
2493:     int a = ptr[Ac];
2494:     if (a == -1) {
2495:       throw std::runtime_error("The node number " + std::to_string(Ac) + 
2496:             " from the spatial values file '" + file_name + " does not belong to the face '" + lFa.name + "'."); 
2497:     }     
2498: 
2499:     lBc.gx[a] = rtmp;
2500:   }
2501: }
2502: 
2503: //-----------------------
2504: // read_temp_spat_values 
2505: //-----------------------
2506: // Read in a file containing temporal and spatial values
2507: // used for a boundary condition.
2508: //
```

```cpp
2797:   gFa.name = "face from traction vtp";
2798:   face_match(com_mod, lFa, gFa, ptr);
2799: 
2800:   // Copy pressure/traction data to MB data structure
2801:   //
2802:   if (lMB.dof == 1) {
2803:     for (int a = 0; a < lFa.nNo; a++) {
2804:       int Ac = ptr[a];
2805:       lMB.d(0,a,0) = -tmpX1[Ac];
2806:       lMB.d(0,a,1) = -tmpX1[Ac];
2807:     }
2808: 
2809:   } else { 
2810:     for (int a = 0; a < lFa.nNo; a++) {
2811:       int Ac = ptr[a];
2812:       for (int i = 0; i < lMB.dof; i++) {
2813:         lMB.d(i,a,0) = tmpX2(i,Ac);
2814:         lMB.d(i,a,1) = tmpX2(i,Ac);
2815:       }
2816:     }
2817:   }
2818: }
2819: 
2820: //-----------
```

## external/flow_solver_source/Code/Source/solver/set_bc.cpp

SHA256: bc1efd7c1a9b3ebc642a6e50fd368567b3f822863068547d6bc8e3082ebbc28b

```cpp
1990:   Array3<double> lK(dof*dof,eNoN,eNoN);
1991: 
1992:   // Constructing LHS/RHS contribution and assembiling them
1993:   //
1994:   for (int e = 0; e < lFa.nEl; e++) {
1995:     cDmn = all_fun::domain(com_mod, com_mod.msh[iM], cEq, lFa.gE(e));
1996:     auto cPhys = eq.dmn[cDmn].phys;
1997: 
1998:     if (lFa.eType == ElementType::NRB) {
1999:       //CALL NRBNNXB(msh(iM), lFa, e)
2000:     }
2001: 
2002:     for (int a = 0; a < eNoN; a++) {
2003:       int Ac = lFa.IEN(a,e);
2004:       ptr(a) = Ac;
2005:       hl.set_col(a, hg.col(Ac));
2006:     }
2007: 
2008:     lK = 0.0;
2009:     lR = 0.0;
2010: 
2011:     for (int g = 0; g < lFa.nG; g++) {
2012:       auto Nx = lFa.Nx.slice(g);
2013:       const Vector<double> nV = nn::gnnb(com_mod, lFa, e, g, Nx, solutions);
2014:       double Jac = utils::norm(nV);
2015:       double w = lFa.w(g)*Jac;
2016:       N = lFa.N.col(g);
2017: 
2018:       Vector<double> h(nsd);
2019:       for (int a = 0; a < eNoN; a++) {
2020:         h = h + N(a)*hl.col(a);
2021:       }
2022: 
2023:       for (int a = 0; a < eNoN; a++) {
2024:         for (int i = 0; i < nsd; i++) {
2025:           lR(i,a) = lR(i,a) - w*N(a)*h(i);
2026:         }
2027:       }
2028:     }
2029: 
2030:     eq.linear_algebra->assemble(com_mod, eNoN, ptr, lK, lR);
2031:   }
2032: }
2033: 
2034: /// @brief Treat Neumann boundaries that are not deforming.
2035: ///
2036: /// Leave the row corresponding to the master node of the owner
2037: /// process in the LHS matrix and the residual vector untouched. For
2038: /// all the other nodes of the face, set the residual to be 0 for
```

## external/flow_solver_source/Code/Source/solver/baf_ini.cpp

SHA256: e4b8946bb52fee5593429490622bc02abd1eac7a1120c6fd5164a3ab52d90680

```cpp
431:   } else if (btest(lBc.bType, enum_int(BoundaryConditionType::bType_ud))) { 
432:     for (int a = 0; a < lFa.nNo; a++) {
433:       int Ac = lFa.gN(a);
434:       s(Ac) = lBc.gx(a);
435:     }
436:   }
437: 
438:   // Now correcting the inlet BC for the inlet ring
439:   //
440:   if (btest(lBc.bType, enum_int(BoundaryConditionType::bType_zp))) { 
441:     for (int jFa = 0; jFa < com_mod.msh[iM].nFa; jFa++) {
442:       if (jFa == iFa) {
443:         continue; 
444:       }
445:       for (int a = 0; a < com_mod.msh[iM].fa[jFa].nNo; a++) {
446:         int Ac = com_mod.msh[iM].fa[jFa].gN(a);
447:         s(Ac) = 0.0;
448:       }
449:     }
450:   }
451: 
452:   // Normalizing the profile for flux
453:   //
454:   double tmp = 1.0;
455:   if (btest(lBc.bType, enum_int(BoundaryConditionType::bType_flx))) {
456:     tmp = all_fun::integ(com_mod, cm_mod, lFa, s, solutions, false);
457:     if (is_zero(tmp)) {
458:       tmp = 1.0;
459:       throw std::runtime_error("Face '" + lFa.name + "' used for a BC has no non-zero node.");
460:      }
461:   }
462: 
463:   for (int a = 0; a < lFa.nNo; a++) {
464:     int Ac = lFa.gN(a);
465:     lBc.gx(a) = s(Ac) / tmp;
466:   }
467: }
468: 
469: //----------
470: // face_ini
471: //----------
472: //
473: void face_ini(Simulation* simulation, mshType& lM, faceType& lFa, const SolutionStates& solutions)
474: {
475:   const auto& Do = solutions.old.get_displacement();
476: 
```
