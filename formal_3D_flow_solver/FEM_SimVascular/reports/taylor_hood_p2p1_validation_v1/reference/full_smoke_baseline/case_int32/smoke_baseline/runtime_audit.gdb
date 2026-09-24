set pagination off
set confirm off
set breakpoint pending on
break *'fs::get_thood_fs(ComMod&, std::array<fsType, 2ul>&, mshType const&, bool, int)'
commands
silent
set $mesh = (char*)$rdx
set $fsp = *(char**)($mesh+896)
printf "TH_RUNTIME mesh_eNoN=%d nFs=%d mesh_eType=%d velocity_eNoN=%d pressure_eNoN=%d velocity_eType=%d pressure_eType=%d vmsStab=%d option=%d\n", *(int*)($mesh+40), *(int*)($mesh+64), *(int*)($mesh+36), *(int*)($fsp+8), *(int*)($fsp+208), *(int*)($fsp+4), *(int*)($fsp+204), (int)$rcx, (int)$r8
disable 1
continue
end
run
