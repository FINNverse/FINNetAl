The zip file Kupers_et_al is part of the manuscript titled:
"Dry season soil water potential maps of a 50 hectare tropical forest plot on Barro Colorado Island, Panama"
The manuscript is published as a Data Descriptor in the journal Scientific Data.
The data in the zip file are accessible through the link in the Data Descriptor.

The input folder in Kupers_et_al includes four files:

1. BCI_Soil_moisture_R_Code.R
	Includes the code to generate soil water potential maps in R.
	
2. BCI_Soil_moisture_mapping.txt
	Includes the following variables:
		id				Unique ID of observation
		location.id		ID of the sampling location
		x				x coordinate of the 50-ha plot
		y				y coordinate of the 50-ha plot
		period			Sampling period (Feb, Mar, Apr 2015 or Mar 2016)
		col.date		Numerical date (i.e. date since 1 January 1900)
		date			Date (M/DD/YYYY)
		col.time		Collection time as fraction of 24 hours
		time			Collection time (hh:mm)
		depth.cat		Depth category. Shallow = 15 cm, middle  = 40 cm, deep = 100 cm
		depth			Sampling depth (cm)
		swp				Soil water potential (MPa)
		swc				Gravimetric soil water content calculated from fresh mass (f) and dry mass (d) as SWC = (f-d)/d 
			
3. BCI_Soil_moisture small_scale.txt
	Includes the following variables:
		id				Unique ID of observation
		location.id		ID of the sampling location
		x				x coordinate of the 50-ha plot
		y				y coordinate of the 50-ha plot
		col.date		Numerical date (i.e. date since 1 January 1900)
		date			Date (M/DD/YYYY)
		col.time		Collection time as fraction of 24 hours
		time			Collection time (hh:mm)
		distance		Distance from the center of the sampling location
		direction		Compass direction relative to the center of the sampling location
		depth			Sampling depth (cm)
		swp				Soil water potential (MPa)
		swc				Gravimetric soil water content calculated from fresh mass (f) and dry mass (d) as SWC = (f-d)/d 

4. BCI_SWP_RF_model.RData
	Includes the rfsrc file for R with the Random Forest model to create soil water potential maps for custom dates.
	
5. BCI_soil_type.txt
	Includes soil type data as digitized from Baillie et al. 2007, emi-detailed soil survey of Barro Colorado Island, Panama. Figure 7.1. https://biogeodb.stri.si.edu/bioinformatics/bci_soil_map/