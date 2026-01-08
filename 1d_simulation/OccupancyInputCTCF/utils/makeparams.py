def paramdict_to_filename(paramdict):

    paramdict_keys={
                    'CTCF_facestall':'face',
                    'CTCF_backstall':'back',
                    'CTCF_lifetime':'Clife',
                    'CTCF_offtime':'Cof',
                    'LEF_lifetime':'life',
                    'LEF_stalled_lifetime':'slife',
                    'LEF_birth':'birth',
                    'LEF_pause':'pause',
                    'LEF_separation':'sep',
                    'sites_per_monomer':'site',
                    'monomers_per_replica':'monomer',
                    'number_of_replica':'replica',
                    'steps':'steps',
                    'velocity_multiplier':'vel'
                    }
    
    filename='file'
    for i in range(len(paramdict)):
        filename += ('_'+paramdict_keys[list(paramdict)[i][:]]+'_'+str(paramdict[list(paramdict)[i]]))
    chars = ['[',']']
    filename_new = ''.join(i for i in filename if not i in chars)    
    return filename_new