 function [constellations] = initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag)
            
            % SYNTAX:
            %   [constellations] = initConstellation(GPS_flag, GLO_flag, GAL_flag, BDS_flag, QZS_flag, SBS_flag);
            %
            % INPUT:
            %   GPS_flag = boolean flag for enabling/disabling GPS usage
            %   GLO_flag = boolean flag for enabling/disabling GLONASS usage
            %   GAL_flag = boolean flag for enabling/disabling Galileo usage
            %   BDS_flag = boolean flag for enabling/disabling BeiDou usage
            %   QZS_flag = boolean flag for enabling/disabling QZSS usage
            %   SBS_flag = boolean flag for enabling/disabling SBAS usage (for ranging)
            %
            % OUTPUT:
            %   constellations = struct with multi-constellation settings
            %
            % DESCRIPTION:
            %   Multi-constellation settings and initialization.
            
            GPS_PRN = [1:32];
            GLO_PRN = [1:24];
            GAL_PRN = [1:30];
            BDS_PRN = [1:37];
            QZS_PRN = [193:196];
            SBS_PRN = 0; %SBAS ranging not supported yet
            
            constellations.GPS     = struct('numSat', numel(GPS_PRN), 'enabled', GPS_flag, 'indexes', 0, 'PRN', GPS_PRN, 'sysID', 'G');
            constellations.GLONASS = struct('numSat', numel(GLO_PRN), 'enabled', GLO_flag, 'indexes', 0, 'PRN', GLO_PRN, 'sysID', 'R');
            constellations.Galileo = struct('numSat', numel(GAL_PRN), 'enabled', GAL_flag, 'indexes', 0, 'PRN', GAL_PRN, 'sysID', 'E');
            constellations.BeiDou  = struct('numSat', numel(BDS_PRN), 'enabled', BDS_flag, 'indexes', 0, 'PRN', BDS_PRN, 'sysID', 'C');
            constellations.QZSS    = struct('numSat', numel(QZS_PRN), 'enabled', QZS_flag, 'indexes', 0, 'PRN', QZS_PRN, 'sysID', 'J');
            constellations.SBAS    = struct('numSat', numel(SBS_PRN), 'enabled', 0,        'indexes', 0, 'PRN', SBS_PRN, 'sysID', 'S'); %SBAS ranging not supported yet
            
            nSatTot = 0; %total number of satellites used given the enabled constellations
            q = 0;       %counter for enabled constellations
            
            systems = fieldnames(constellations);
            constellations.indexes = [];
            constellations.PRN = [];
            constellations.systems = [];
            for i = 1 : numel(systems)
                if(constellations.(systems{i}).enabled)
                    nSatTot = nSatTot + constellations.(systems{i}).numSat;
                    q = q + 1;
                    if (q == 1)
                        indexes_tmp = [1 : constellations.(systems{i}).numSat];
                    else
                        indexes_tmp = [indexes_tmp(end) + 1 : indexes_tmp(end) + constellations.(systems{i}).numSat];
                    end
                    constellations.(systems{i}).indexes = indexes_tmp;
                    constellations.indexes = [constellations.indexes, indexes_tmp];
                    constellations.PRN = [constellations.PRN, constellations.(systems{i}).PRN];
                    constellations.systems = [constellations.systems char(ones(size(constellations.(systems{i}).indexes))*constellations.(systems{i}).sysID)];
                end
            end
            
            constellations.nEnabledSat = nSatTot;
        end
        