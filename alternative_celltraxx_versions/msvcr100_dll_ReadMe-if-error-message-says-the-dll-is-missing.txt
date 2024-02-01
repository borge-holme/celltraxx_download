In some cases, celltraxx.exe may not run even if the correct Microsoft Visual C++ Redistributable has been installed. Then starting celltraxx.exe may give an error message which claims that the msvcr100.dll file is missing, even if it is present in the C:\Windows\System32\ directory. In such cases that file may be too old. 

It has helped to copy the msvcr100.dll from this folder (which was created 2011-06-11) directly into the C:\celltraxx_system\ directory. The file provided here is the 64-bit version (x64). If using a 32-bit Windows system (x86), please download the corresponding 32-bit version from the net. See for instance:

https://answers.microsoft.com/en-us/windows/forum/all/msvcr100-dll/826a842c-7a2f-4ebc-bb21-4cee86b234b0

