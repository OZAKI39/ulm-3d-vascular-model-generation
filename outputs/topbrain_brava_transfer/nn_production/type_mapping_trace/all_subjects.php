<html lang="en">
<head>

<script type="text/javascript">
// Javascript function *****************************************************************************************************
function metric_function(link)
{
	var metric_name=link[link.selectedIndex].value;

	if (metric_name == '');
	else
		location.href = "all_subjects.php?+name_metric="+metric_name+"&ad_search=";			
}
</script>

<meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">    
<title>Browse all files</title>
<link rel="stylesheet" type="text/css" href="sheet1.css">
</head>
<body>
<div align="center">
<script language="JavaScript">
<!-- This script calculate the odiern date -->
data = new Date();
day = data.getDay();
month = data.getMonth();
date= data.getDate();
year= data.getYear();
if(year<1900)year=year+1900;
if(day == 0) day = " Sunday, ";
if(day == 1) day = " Monday, ";
if(day == 2) day = " Tuesday, ";
if(day == 3) day = " Wednesday, ";
if(day == 4) day = " Thursday, ";
if(day == 5) day = " Friday, ";
if(day == 6) day = " Saturday, ";
if(month == 0) month = "January ";
if(month ==1) month = "February ";
if(month ==2) month = "March ";
if(month ==3) month = "April ";
if(month ==4) month = "May ";
if(month ==5) month = "June ";
if(month ==6) month = "July ";
if(month ==7) month = "August ";
if(month ==8) month = "September ";
if(month ==9) month = "October ";
if(month ==10) month = "November ";
if(month ==11) month = "December";
</script>
<table border="0" width="1000" cellpadding="0" cellspacing="0">
<tr>
   <td width="20%" align="center" bgcolor="#038100">     
 	    <table border="0" width="960" cellpadding="0" cellspacing="0" height="2">  		
		<tr>
           <td width="100%" align="center" bgcolor="#038100"> 
           </td>   
        </tr>
   		</table>
        <table border="0" width="1000" cellpadding="0" cellspacing="0" height="120">
        <tr>
           <td width="100%" align="center"> 
            <img src="images/title1.png" height="110" width="1000">
           </td>   
        </tr>
        </table>
	    <table border="0" width="960" cellpadding="0" cellspacing="0" height="20">
        <tr>
           <td width="30%" align="left" bgcolor="#038100"> 
            <font face="Verdana, Arial, Helvetica, sans-serif" size='2' color="#FFFFFF"><strong>
             <script language="JavaScript">   
               document.write(""+day+" "+month+" "+date+" "+year+"");  
             </script>   
            </strong>
	    </font>
           </td>
           <td width="40%" align="center" bgcolor="#038100">
            	<hr width="90%" color="#FFFFFF">
           </td>
            <td width="30%" align="right" bgcolor="#038100"> 
            <font face="Verdana, Arial, Helvetica, sans-serif" size='2' color="#FFFFFF"><strong>Version 1.1</strong></font>
           </td>  
        </tr>
        </table>
   </td>
</tr>
</table><table width="1000" border="0" cellpadding="0" cellspacing="0">
<tr>
	<td width="200" bgcolor="#FFFFFF" valign="top">
     <br>    
	<br><center><div id='scelta_menu'><a href='index.html' style='text-decoration:none'><font color='#006633' >Home page</font></a><br><hr width='90%'></div><div id='scelta_menu'><a href='all_subjects.php?clear=1' style='text-decoration:none'><font color='#A00000'>Browse all files</font></a><br><hr width='90%'></div><div id='scelta_menu'><a href='search_advanced.php' style='text-decoration:none'><font color='#006633' >Search files</font></a><br><hr width='90%'></div><div id='scelta_menu'><a href='software.php' style='text-decoration:none'><font color='#006633' >Related Links</font></a><br><hr width='90%'></div><div id='scelta_menu'><a href='terms_of_use.php' style='text-decoration:none'><font color='#006633'>Terms of use</font></a><br><hr width='90%'></div><div id='scelta_menu'><a href='help.php' style='text-decoration:none'><font color='#006633' >Help</font></a><br><hr width='90%'></div>	
	</td>
	<td width="800" align="center" valign="top" >
	<!-- BODY -->	
	<div id='div_menu6'>
	<br>
	<font size="2" color="#000099" face="Verdana, Arial, Helvetica, sans-serif">&nbsp;
	</font>
	<BR>
		<table width="85%" border='0'>
		<tr>
			<td width="20%" align="center">
			</td>
			<td width="60%" align="center">
			
			<div id='td1'><font size='3' color='white'>			
			61 subjects are present in the database			</font>
			</td>	
			</div>
			<td width="20%" align="center">
			</td>	
		</tr>
	</table>	

	<br>

	<div id='div_menu1'>

		<table width="85%" border='0' cellspacing="5">
		<tr>
			<td width="20%" align="center">
				<input type='button' value='Download swc files' onClick="location.href='files/swc_files.zip'" id='button1'>			</td>
			<td width="60%" align="center">
			<form action='metrics.php' method="post" style="display:inline" target="_blank">
				<input type="submit" name='metrics' value='View statistics for metrics' id='button1'>
			<input type='hidden' name='id_serial' value='a:61:{i:0;s:1:"1";i:1;s:1:"2";i:2;s:1:"3";i:3;s:1:"4";i:4;s:1:"5";i:5;s:1:"6";i:6;s:1:"7";i:7;s:1:"8";i:8;s:1:"9";i:9;s:2:"10";i:10;s:2:"11";i:11;s:2:"12";i:12;s:2:"13";i:13;s:2:"14";i:14;s:2:"15";i:15;s:2:"16";i:16;s:2:"17";i:17;s:2:"18";i:18;s:2:"20";i:19;s:2:"21";i:20;s:2:"22";i:21;s:2:"23";i:22;s:2:"24";i:23;s:2:"25";i:24;s:2:"26";i:25;s:2:"27";i:26;s:2:"28";i:27;s:2:"29";i:28;s:2:"30";i:29;s:2:"31";i:30;s:2:"32";i:31;s:2:"33";i:32;s:2:"34";i:33;s:2:"35";i:34;s:2:"36";i:35;s:2:"37";i:36;s:2:"38";i:37;s:2:"39";i:38;s:2:"40";i:39;s:2:"41";i:40;s:2:"42";i:41;s:2:"43";i:42;s:2:"44";i:43;s:2:"45";i:44;s:2:"46";i:45;s:2:"47";i:46;s:2:"48";i:47;s:2:"49";i:48;s:2:"50";i:49;s:2:"51";i:50;s:2:"52";i:51;s:2:"53";i:52;s:2:"54";i:53;s:2:"55";i:54;s:2:"56";i:55;s:2:"57";i:56;s:2:"58";i:57;s:2:"59";i:58;s:2:"60";i:59;s:2:"61";i:60;s:2:"62";}'>
			</form>			
			</td>	
			<td width="20%" align="center">		
			</td>	
		</tr>
		</table>
	
	
		<table width="80%" border='0' cellpadding="0" cellspacing="2">
		<tr>
			<td width="75%" align="center" id='td2'>
			<input type='radio' name='type_metric' value='' onClick="javascript:location.href='all_subjects.php?type_metric=complete&ad_search='" checked='checked'>Complete<input type='radio' name='type_metric' value='' onClick="javascript:location.href='all_subjects.php?type_metric=aca&ad_search='">ACA<input type='radio' name='type_metric' value='' onClick="javascript:location.href='all_subjects.php?type_metric=mca&ad_search='">MCA<input type='radio' name='type_metric' value='' onClick="javascript:location.href='all_subjects.php?type_metric=pca&ad_search='">PCA		
			</td>	
			<td width="5%" align="center">
			</td>
			<td width="20%" align="center">
				<select name='metric' size='1' cols='1' onChange="metric_function(this)">
				<option value=''>Choose a metric</option>	
				<option value=''>-------</option>						
				<option value='total_n_branches'>Total number of branches</option>
				<option value='total_lenght'>Total length</option>
				<option value='max_branch_order'>Max branch order</option>
				<option value='max_path_distance'>Max path distance</option>
				<option value='max_euclidean_distance'>Max euclidean distance</option>
				<option value='av_width'>Width</option>
				<option value='av_height'>Height</option>
				<option value='av_depth'>Depth</option>
				<option value='av_bifurcation_amplitude_local'>Mean local bifurcation amplitude</option>
				<option value='av_bifurcation_amplitude_remote'>Mean bifurcation amplitude</option>
				<option value='av_contraction'>Contraction</option>
				<option value='av_partition_asymmetry'>Partition asymmetry</option>
				<option value='av_fractal_dimension'>Fractal dimension</option>
				<option value='av_branch_path_length'>Branch path length</option>
				<option value='av_bifurcation_tilt_remote'>Mean bifurcation tilt</option>
				<option value='av_bifurcation_torque_remote'>Mean bifurcation torque</option>
				<option value='av_fragmentation'>Mean fragmentation</option>
				</select> 		
			</td>	
		</tr>
		</table>
		<br>
	</div>
	
	<br>
	<div id='table1'>
				<table border='0' width='75%'>
			<tr>
				<td width='60%' align='center'>
				<a href='all_subjects.php?order=id_code'><font color='white' size='2'>CODE</a> </font> </td>
				<td width='20%' align='center'>
				<a href='all_subjects.php?order=age'><font color='white' size='2'> AGE </a></font></td>
				<td width='20%' align='center'>
				<a href='all_subjects.php?order=sex'><font color='white' size='2'> SEX </a></font></td>     
			</tr>
			</table> 
			  
	</div> 

			<table border='0' width='75%'>
			
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=1' target='_blank'><font size='2' id='font1'>BG001</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 24</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=2' target='_blank'><font size='2' id='font1'>BG0002</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 31</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=3' target='_blank'><font size='2' id='font1'>BG0003</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 29</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=4' target='_blank'><font size='2' id='font1'>BG04</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 21</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=5' target='_blank'><font size='2' id='font1'>BG05</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 35</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=6' target='_blank'><font size='2' id='font1'>BG06</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 20</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=7' target='_blank'><font size='2' id='font1'>BG07</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 41</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=8' target='_blank'><font size='2' id='font1'>BG08</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 27</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=9' target='_blank'><font size='2' id='font1'>BG09</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 23</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=10' target='_blank'><font size='2' id='font1'>BG10</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 21</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=11' target='_blank'><font size='2' id='font1'>BG11</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 44</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=12' target='_blank'><font size='2' id='font1'>BG12</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 36</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=13' target='_blank'><font size='2' id='font1'>BG13</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 38</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=14' target='_blank'><font size='2' id='font1'>BG15</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 27</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=15' target='_blank'><font size='2' id='font1'>BG17</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 31</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=16' target='_blank'><font size='2' id='font1'>BG18</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 24</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=17' target='_blank'><font size='2' id='font1'>Set 8</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 32</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=18' target='_blank'><font size='2' id='font1'>Set 9</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 46</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=20' target='_blank'><font size='2' id='font1'>BG0014</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 23</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=21' target='_blank'><font size='2' id='font1'>BG0019</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 34</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=22' target='_blank'><font size='2' id='font1'>BG0020</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 31</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=23' target='_blank'><font size='2' id='font1'>BG0021</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 33</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=24' target='_blank'><font size='2' id='font1'>BG0022</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 21</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=25' target='_blank'><font size='2' id='font1'>BH0003</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 64</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=26' target='_blank'><font size='2' id='font1'>BH0004</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 19</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=27' target='_blank'><font size='2' id='font1'>BH0005</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 28</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=28' target='_blank'><font size='2' id='font1'>BH0006</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 26</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=29' target='_blank'><font size='2' id='font1'>BH0008</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 28</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=30' target='_blank'><font size='2' id='font1'>BH0009</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 42</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=31' target='_blank'><font size='2' id='font1'>BH0010</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 26</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=32' target='_blank'><font size='2' id='font1'>BH0011</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 27</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=33' target='_blank'><font size='2' id='font1'>BH0012</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 43</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=34' target='_blank'><font size='2' id='font1'>BH0013</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 28</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=35' target='_blank'><font size='2' id='font1'>BH0014</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 23</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=36' target='_blank'><font size='2' id='font1'>BH0015</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 25</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=37' target='_blank'><font size='2' id='font1'>BH0016</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 28</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=38' target='_blank'><font size='2' id='font1'>BH0017</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 24</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=39' target='_blank'><font size='2' id='font1'>BH0018</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 21</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=40' target='_blank'><font size='2' id='font1'>BH0019</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 27</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=41' target='_blank'><font size='2' id='font1'>BH0020</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 21</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=42' target='_blank'><font size='2' id='font1'>BH0021</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 33</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=43' target='_blank'><font size='2' id='font1'>BH0022</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 47</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=44' target='_blank'><font size='2' id='font1'>BH0023</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 19</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=45' target='_blank'><font size='2' id='font1'>BH0024</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 22</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=46' target='_blank'><font size='2' id='font1'>BH0025</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 24</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=47' target='_blank'><font size='2' id='font1'>BH0026</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 37</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=48' target='_blank'><font size='2' id='font1'>BH0027</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 46</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=49' target='_blank'><font size='2' id='font1'>BH0029</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 46</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=50' target='_blank'><font size='2' id='font1'>BH0030</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 33</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=51' target='_blank'><font size='2' id='font1'>BH0031</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 24</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=52' target='_blank'><font size='2' id='font1'>BH0032</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 48</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=53' target='_blank'><font size='2' id='font1'>BH0033</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 59</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=54' target='_blank'><font size='2' id='font1'>BH0034</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 27</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=55' target='_blank'><font size='2' id='font1'>BH0035</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 24</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=56' target='_blank'><font size='2' id='font1'>BH0036</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 29</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=57' target='_blank'><font size='2' id='font1'>BH0037</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 50</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=58' target='_blank'><font size='2' id='font1'>BH0038</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 22</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=59' target='_blank'><font size='2' id='font1'>BH0039</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 42</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=60' target='_blank'><font size='2' id='font1'>BH0040</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 30</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb2'>
							<a href='query_subject.php?id=61' target='_blank'><font size='2' id='font1'>BI0001</a></font></td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> 22</font> </td>
							<td width='20%' align='center' id='tb2'><font color='black' size='2'> M</font> </td>
							</tr>
							
							<tr>
							<td width='60%' align='center' id='tb5'>
							<a href='query_subject.php?id=62' target='_blank'><font size='2' id='font1'>BH0002</a></font></td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> 25</font> </td>
							<td width='20%' align='center' id='tb5'><font color='black' size='2'> F</font> </td>
							</tr>
							       
		</table>
	
	
	<br><br><br>
	</td>
		
	</div>
</tr>
</table>

</div>
</body>
</html>