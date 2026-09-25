<html lang="en">
<head>
<title>Help</title>
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
</table><table width="1000" height="400" border="0" cellpadding="0" cellspacing="0">
<tr>
	<td width="200" bgcolor="#FFFFFF" valign="top">
        <br>    
	<br><center><div id='scelta_menu'><a href='index.html' style='text-decoration:none'><font color='#006633' >Home page</font></a><br><hr width='90%'></div><div id='scelta_menu'><a href='all_subjects.php?clear=1' style='text-decoration:none'><font color='#006633' >Browse all files</font></a><br><hr width='90%'></div><div id='scelta_menu'><a href='search_advanced.php' style='text-decoration:none'><font color='#006633' >Search files</font></a><br><hr width='90%'></div><div id='scelta_menu'><a href='software.php' style='text-decoration:none'><font color='#006633' >Related Links</font></a><br><hr width='90%'></div><div id='scelta_menu'><a href='terms_of_use.php' style='text-decoration:none'><font color='#006633'>Terms of use</font></a><br><hr width='90%'></div><div id='scelta_menu'><a href='help.php' style='text-decoration:none'><font color='#A00000'>Help</font></a><br><hr width='90%'></div>	</td>
	<td width="800" align="center" valign="top">
	<!-- BODY -->
	<div id='div_menu7'>
		<br><br>
		<table border="0" width="40%">
		<tr>
			<td width='100%' align="center">
			
							<div id='td1'>
								<font size='4'>Help</font>
							</div>		
							<br>					
								
			</td>
		</tr>
		</table>
		<br>
		<table width="85%" border="0">
		<tr>
			<td width="100%">
			<font face="Verdana, Arial, Helvetica, sans-serif" size='3'>
				
			<font color="#006200">
			The 'Browse' tab displays all arterial reconstructions, linking each to its individual record. Every record gives access to the digital reconstruction (swc) file and includes  a 3D animation and a web browser-based interactive display. The 'Search' tab allows users to find records based on the subject, gender, or age, and by the reconstruction morphometrics.The de-identified MRA image stacks for all records are available for <a href='http://www.mdbrain.org/personalpages/peter/mra_brains.tar' target='_blank'>bulk download</a>.<br><br>
			Please address questions or comments to <a href='mailto:ascoli@gmu.edu' id='email_ascoli'><strong>ascoli@gmu.edu</strong></a>
			</font>
	</div>	
	</td>
	</tr>
	</table>
	
    </td>
</tr>
</table>
</div>
</body>
</html>
