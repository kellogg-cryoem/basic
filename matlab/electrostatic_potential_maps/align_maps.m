function[aligned] = align_maps(test,ref)
    thresh = 0.001;
    test_ = test;
    test_(find(test_ < thresh)) = 0;
    ref_ = ref;
    ref_(find(ref_ < thresh)) = 0;
    testcom = com_map(test_);
    refcom = com_map(ref_);
    refcom-testcom
    aligned = imtranslate(test,(refcom-testcom));

end