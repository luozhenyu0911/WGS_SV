def combine_line(f, baseline, pass_line, cutoff = 0.5):
    try:
        line =  next(f)
    except:
        pass_line.append(baseline.strip())
        return pass_line
    END, SVLEN, CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT, *gt_list = line_parser(baseline)
    END2, SVLEN2, CHROM2,POS2,ID2,REF2,ALT2,QUAL2,FILTER2,INFO2,FORMAT2, *gt_list2 = line_parser(line)
    if POS2 >= END or CHROM2 != CHROM or ALT2 != ALT:
        pass_line.append(baseline.strip())
        baseline = line
    elif POS <= POS2 < END:
        if END2 <= END:
            overlap_length = END2 - POS2
            
        elif END <= END2:
            overlap_length = END - POS2
            
        if (cutoff*int(SVLEN2) <= overlap_length) or (cutoff*int(SVLEN) <= overlap_length):
            maxEND = max(END, END2)
            maxSVLEN = max(END, END2) - POS
            new_gt_list = combine_gt(gt_list, SVLEN, gt_list2, SVLEN2)
            new_info = "END=" + str(maxEND) + ";SVLEN=" + str(maxSVLEN) + ";SVTYPE=" + ALT
            new_line = '\t'.join([CHROM, str(POS), ID, REF, ALT, QUAL, FILTER, new_info, FORMAT, '\t'.join(new_gt_list)])
            # pass_line.append(new_line)
            baseline = new_line
        else:
            pass_line.append(baseline.strip())
            baseline = line
    return combine_line(f, baseline, pass_line, cutoff = 0.5)