import configparser

# This is used by the dinara scripts to access the config file.

# Read the config file.
def getConfig():
    config = configparser.ConfigParser()
    configFileName = 'dinara.conf'
    if not config.read(configFileName):
        raise Exception('Error reading config file %s.' % configFileName)
    return config


