# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.
def run():
    import gpbo.core.optimize

    #import gpbo.configs.randomsh as optconfig
    #import gpbo.configs.gridsh as optconfig
    import gpbo.configs.EIMLsh as optconfig
    #import gpbo.configs.PESbssh as optconfig
    #import gpbo.configs.EIpumppid as optconfig
    #import gpbo.configs.PESpump as optconfig
    #import gpbo.configs.PBpump as optconfig

    O = gpbo.search(optconfig)
    print "RESULT: {}".format(O)
    return

if __name__=="__main__":
    run()