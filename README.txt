This project demonstrates how to compute the Allan Variance from time series data.

Files:
	AllanVar_Example.ipynb	: A jupyter notebook containing an example calculation.
	tdmsLib_vObj.py 	: The above jupyter notebook uses functions from this library for the computation.
	tser_class.py 		: Defines the "TSer" class for working with time series data.
	requirements.txt 	: The project should run in a virtual environment that contains these Python libraries.
					To install:
						1) Create a virtual environment for the project.
						2) Enter the environment
						3) Navigate to directory containing "requirements.py".
						4) Install with the command: pip install -r requirements.txt
	
	The raw data for the example is stored in the "ExampleData" folder.
		Two data sets are provided: "GapData_Janike" shows less stable data; "StableChristmasData_Sergii" is more stable.
		Note: The two data sets were also collected with 0.01ms and 0.001ms time steps respectively.
