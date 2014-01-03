ieeg-2edf
=========

Methods to stream data from the IEEG-Portal to an EDF file. We encourage IEEG-Portal users to use the IEEG-Toolbox to access data as this is faster and more reliable. However, it is possible to download data from the portal and stream it to an EDF file using these methods.

This method utilizes the IEEG-Portal Matlab toolbox to connect to the IEEG-Portal and streams the contents of a dataset to a local EDF file which can be opened by 3rd party EEG software packages. 

       Example: Stream first 100000 values of all channels in 'Study 005' to EDF.
           >> session = IEEGSession('Study 005', 'username', 'pwdFile');
           >> ieeg2edf(session.data, [1 100000]);   
