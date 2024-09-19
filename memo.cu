// Copy function for SigArr
    copy3DarrayHostToDevice(befaft_host->sa.Txx , &(befaft_device.sa.Txx) , ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(befaft_host->sa.Txxx, &(befaft_device.sa.Txxx), ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(befaft_host->sa.Txxy, &(befaft_device.sa.Txxy), ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(befaft_host->sa.Txxz, &(befaft_device.sa.Txxz), ran.sr.Txx.x, ran.sr.Txx.y, ran.sr.Txx.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tyy , &(befaft_device.sa.Tyy) , ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tyyx, &(befaft_device.sa.Tyyx), ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tyyy, &(befaft_device.sa.Tyyy), ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tyyz, &(befaft_device.sa.Tyyz), ran.sr.Tyy.x, ran.sr.Tyy.y, ran.sr.Tyy.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tzz , &(befaft_device.sa.Tzz) , ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tzzx, &(befaft_device.sa.Tzzx), ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tzzy, &(befaft_device.sa.Tzzy), ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);
    copy3DarrayHostToDevice(befaft_host->sa.Tzzz, &(befaft_device.sa.Tzzz), ran.sr.Tzz.x, ran.sr.Tzz.y, ran.sr.Tzz.z);

    // Copy function for TauArr
    copy3DarrayHostToDevice(befaft_host->ta.Txy , &(befaft_device.ta.Txy) , ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    copy3DarrayHostToDevice(befaft_host->ta.Txyx, &(befaft_device.ta.Txyx), ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    copy3DarrayHostToDevice(befaft_host->ta.Txyy, &(befaft_device.ta.Txyy), ran.tr.Txy.x, ran.tr.Txy.y, ran.tr.Txy.z);
    copy3DarrayHostToDevice(befaft_host->ta.Tyz , &(befaft_device.ta.Tyz) , ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    copy3DarrayHostToDevice(befaft_host->ta.Tyzy, &(befaft_device.ta.Tyzy), ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    copy3DarrayHostToDevice(befaft_host->ta.Tyzz, &(befaft_device.ta.Tyzz), ran.tr.Tyz.x, ran.tr.Tyz.y, ran.tr.Tyz.z);
    copy3DarrayHostToDevice(befaft_host->ta.Tzx , &(befaft_device.ta.Tzx) , ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);
    copy3DarrayHostToDevice(befaft_host->ta.Tzxz, &(befaft_device.ta.Tzxz), ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);
    copy3DarrayHostToDevice(befaft_host->ta.Tzxx, &(befaft_device.ta.Tzxx), ran.tr.Tzx.x, ran.tr.Tzx.y, ran.tr.Tzx.z);

    // Copy function for VelArr
    copy3DarrayHostToDevice(befaft_host->va.Vx , &(befaft_device.va.Vx) , ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    copy3DarrayHostToDevice(befaft_host->va.Vxx, &(befaft_device.va.Vxx), ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    copy3DarrayHostToDevice(befaft_host->va.Vxy, &(befaft_device.va.Vxy), ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    copy3DarrayHostToDevice(befaft_host->va.Vxz, &(befaft_device.va.Vxz), ran.vr.Vx.x, ran.vr.Vx.y, ran.vr.Vx.z);
    copy3DarrayHostToDevice(befaft_host->va.Vy , &(befaft_device.va.Vy) , ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    copy3DarrayHostToDevice(befaft_host->va.Vyx, &(befaft_device.va.Vyx), ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    copy3DarrayHostToDevice(befaft_host->va.Vyy, &(befaft_device.va.Vyy), ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    copy3DarrayHostToDevice(befaft_host->va.Vyz, &(befaft_device.va.Vyz), ran.vr.Vy.x, ran.vr.Vy.y, ran.vr.Vy.z);
    copy3DarrayHostToDevice(befaft_host->va.Vz , &(befaft_device.va.Vz) , ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    copy3DarrayHostToDevice(befaft_host->va.Vzx, &(befaft_device.va.Vzx), ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    copy3DarrayHostToDevice(befaft_host->va.Vzy, &(befaft_device.va.Vzy), ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);
    copy3DarrayHostToDevice(befaft_host->va.Vzz, &(befaft_device.va.Vzz), ran.vr.Vz.x, ran.vr.Vz.y, ran.vr.Vz.z);