from HinetPy import Client, win32
import os

# You need a Hi-net account to access their data
# User and password must be saved as enviroment variables in the .bashrc file (linux)
# Example
# export USER_HiNet='juanito'
# export PASSWORD_HiNet='secreto' 

client = Client(os.environ['USER_HiNet'], os.environ['PASSWORD_HiNet'])
client.select_stations("0101", ["N.KWSH", "N.IWNH","N.ASAH","N.NKWH", "N.YHTH"])
data, ctable = client.get_continuous_waveform('0101', '201001010000', 14400)
# The request and downloading process will take several minutes
# waiting data request ...
# waiting data downloading ...

win32.extract_sac(data, ctable)


win32.extract_sacpz(ctable)

