# UI Debug: LabFlow å‰ç«¯

## è®¾è®¡å†³ç­–è®°å½•

### [2026-02-11] ç‚¹å‡»æŒ‰é’®åçš„å³æ—¶åé¦ˆä¸è¿è¡Œæ—¥å¿—ã€Œå°ç”µè§†ã€

- **ç”¨æˆ·éœ€æ±‚**ï¼šæ¯ä¸€æ­¥ç‚¹æŒ‰é’®åè¦æœ‰å³æ—¶åé¦ˆï¼›è®­ç»ƒç­‰é•¿ä»»åŠ¡è¦æœ‰å®æ—¶ç›‘æµ‹ï¼Œå¦‚è¿›åº¦æ¡æˆ–æŠŠåç«¯ shell å½“ã€Œå°ç”µè§†ã€å±•ç¤ºã€‚
- **æ–¹æ¡ˆ**ï¼š
  1. **å…¨å±€ã€Œå½“å‰ä»»åŠ¡ã€åé¦ˆ**ï¼šå½“å­˜åœ¨è¿›è¡Œä¸­æ­¥éª¤ï¼ˆ`busyStep`ï¼‰æˆ–ä»»åŠ¡å¤„äº queued/running æ—¶ï¼Œåœ¨ç³»ç»ŸçŠ¶æ€ä¸‹æ–¹å±•ç¤ºä¸€ä¸ªå›ºå®šæ¡å¸¦ï¼ŒåŒ…å«ï¼š
     - æ— é™å¾ªç¯åŠ¨ç”»çš„è¿›åº¦æ¡ï¼ˆindeterminateï¼‰ï¼Œè¡¨ç¤ºã€Œæ­£åœ¨å¤„ç†ã€ï¼›
     - æ–‡æ¡ˆï¼šä¸Šä¼ ä¸­â€¦ / æ ¡éªŒä¸­â€¦ / è§£æ Seurat ä¸­â€¦ / 500Ã—500 é¢„å¤„ç†ä¸­â€¦ / æäº¤è®­ç»ƒä¸­â€¦ / è®­ç»ƒè¿è¡Œä¸­â€¦ã€‚
  2. **è¿è¡Œæ—¥å¿—ã€Œå°ç”µè§†ã€**ï¼šåœ¨ã€Œ6) ä»»åŠ¡è½®è¯¢çŠ¶æ€ã€ä¸­ï¼Œå½“å­˜åœ¨ job æ—¶å±•ç¤ºã€Œè¿è¡Œæ—¥å¿—ã€åŒºå—ï¼š
     - è½®è¯¢ `GET /api/jobs/{job_id}/log`ï¼Œæ¯ 2.5 ç§’åˆ·æ–°ï¼›
     - å†…å®¹åœ¨æ·±è‰²ã€ç­‰å®½ã€å¯æ»šåŠ¨çš„ `<pre>` ä¸­å±•ç¤ºï¼ˆmax-height 220pxï¼‰ï¼Œé£æ ¼ç±»ä¼¼ç»ˆç«¯ï¼›
     - ä»»åŠ¡ç»“æŸåä¿ç•™æœ€åä¸€æ¬¡æ‹‰å–çš„æ—¥å¿—ï¼Œä¸ä¸‹æ–¹ã€ŒæŸ¥çœ‹æ—¥å¿—ã€ä¸€è‡´ã€‚
- **ä»£ç **ï¼š
  - `frontend/src/App.tsx`ï¼šæ–°å¢ `liveJobLog` çŠ¶æ€ï¼›è½®è¯¢ `getJobLog` çš„ useEffectï¼ˆjob ä¸º queued/running æ—¶ï¼‰ï¼›`taskLabelMap`ã€`showTaskFeedback`ã€`taskLabel`ï¼›æ’å…¥ã€Œå½“å‰ä»»åŠ¡ã€æ¡å¸¦ä¸ã€Œè¿è¡Œæ—¥å¿—ã€é¢æ¿ã€‚
  - `frontend/src/styles/tokens.css`ï¼š`.task-feedback`ã€`.task-progress`ï¼ˆå« `task-progress-shift` åŠ¨ç”»ï¼‰ã€`.task-label`ã€`.live-log-box`ã€`.live-log-pre`ã€‚
- **Checkfix**ï¼š`npm run lint`ã€`npm run build` é€šè¿‡ã€‚
- **åç»­å¯åš**ï¼šè‹¥åç«¯ä¸º validate/inspect/prepare æä¾›è¿›åº¦æˆ–æµå¼æ—¥å¿—ï¼Œå¯å†æ¥ä¸Šè¿›åº¦æˆ–å®æ—¶æ—¥å¿—ã€‚

### [2026-02-11] æ ¡éªŒé”™è¯¯å°±è¿‘å±•ç¤º + R/conda å‚æ•°èœå•åŒ–
- **ç”¨æˆ·éœ€æ±‚**ï¼šæŠ¥é”™åº”åœ¨ç”¨æˆ·è§†çº¿å¤„ï¼ˆæŒ‰é’®æ—ï¼‰å¼¹å‡ºå¹¶ç»™å‡ºå½“å‰ç³»ç»Ÿæ¨èå‚æ•°ï¼›R ç›¸å…³å‚æ•°ï¼ˆconda.batã€R ç¯å¢ƒåï¼‰åº”ç”±åç«¯æ¢æµ‹ååœ¨å‰ç«¯ç”¨èœå•é€‰æ‹©ï¼Œé¿å…æ‰‹è¾“ã€‚
- **å®ç°**ï¼šâ‘  æ ¡éªŒå¤±è´¥æ—¶ä»…è®¾ç½® validateErrorï¼ˆä¸è®¾ globalErrorï¼‰ï¼Œåœ¨ã€Œå¼€å§‹æ ¡éªŒã€æŒ‰é’®æ—ç”¨ step-error åŒºå—å±•ç¤ºé”™è¯¯ä¸ validateRecommendationï¼ˆRscript/Conda æ—¶æ¨è cmd_conda + å®Œæ•´è·¯å¾„ + r-4.3ï¼‰ï¼›è¡¨å•é¡¹ onChange æ—¶æ¸…é™¤ validateErrorã€‚â‘¡ åç«¯æ–°å¢ GET /api/runtime/conda-envsï¼ˆbackend/app/api/runtime.pyï¼‰ï¼šWindows ä¸‹æ‰«æ PATH ä¸å¸¸è§è·¯å¾„å¾—åˆ° conda.bat å€™é€‰ï¼Œæ‰§è¡Œ conda env listï¼ˆ--json æˆ–è§£ææ–‡æœ¬ï¼‰å¾—åˆ°ç¯å¢ƒååˆ—è¡¨ï¼›æ”¯æŒ query conda_bat åªè¿”å›è¯¥è·¯å¾„ä¸‹çš„ envsã€‚â‘¢ å‰ç«¯ï¼šè¯·æ±‚ getCondaEnvs() å¡«å…… condaBatCandidatesã€condaEnvsListï¼›conda.bat ç”¨ selectï¼ˆå€™é€‰ +ã€Œå…¶ä»–ï¼ˆæ‰‹åŠ¨è¾“å…¥ï¼‰ã€ï¼‰+ å¯é€‰ inputï¼›Conda R ç¯å¢ƒåç”¨ selectï¼ˆconda_envsï¼‰ï¼›åˆ‡æ¢ conda.bat æ—¶é‡æ–°è¯·æ±‚ envsã€‚â‘£ æ ·å¼ï¼š.step-actionsã€.step-errorã€.step-error-msgã€.step-error-recommendã€‚
- **Checkfix**ï¼šruffã€frontend lint/build é€šè¿‡ã€‚

### [2026-02-11] å‚æ•°ååã€Œ?ã€æ‚¬åœ 1 ç§’æ˜¾ç¤ºè¯´æ˜æ°”æ³¡ï¼ˆå‚»ç“œåŒ–ï¼‰
- **ç”¨æˆ·éœ€æ±‚**ï¼šæ¯ä¸ªå‚æ•°åç§°åé¢æ”¾ä¸€ä¸ªå°é—®å·ï¼Œé¼ æ ‡æ‚¬åœçº¦ 1 ç§’åå¼¹å‡ºæ°”æ³¡ï¼Œè§£é‡Šè¯¥å‚æ•°åœ¨å•ç»†èƒæ•°æ®é‡Œæ˜¯ä»€ä¹ˆã€å¹²ä»€ä¹ˆç”¨ï¼ŒæŠŠç”¨æˆ·å½“éä¸“ä¸šç”¨æˆ·ã€Œå–‚åˆ°å˜´é‡Œã€ã€‚
- **æ–¹æ¡ˆ**ï¼š
  1. **ParamTooltip ç»„ä»¶**ï¼ˆ`frontend/src/components/ParamTooltip.tsx`ï¼‰ï¼šæ¸²æŸ“å†…è”ã€Œ?ã€å›¾æ ‡ï¼›`onMouseEnter` è®¾ 1000ms å®šæ—¶å™¨ï¼Œ`onMouseLeave` æ¸…é™¤å¹¶éšè—ï¼›å®šæ—¶åˆ°åæ˜¾ç¤ºæ°”æ³¡ï¼ˆ`role="tooltip"`ï¼‰ï¼Œæ°”æ³¡å†…å¯ç»§ç»­æ‚¬åœä»¥ä¿æŒæ˜¾ç¤ºï¼›æ°”æ³¡æ·±è‰²èƒŒæ™¯ã€ç™½å­—ã€å°ä¸‰è§’æŒ‡å‘è§¦å‘ç‚¹ï¼Œmax-width 320pxã€‚
  2. **æ–‡æ¡ˆ**ï¼šåŒä¸€æ–‡ä»¶å†…å¯¼å‡º `PARAM_TOOLTIPS`ï¼ˆRecord<string, string>ï¼‰ï¼Œä¸ºã€Œæ•°æ®åç§°ã€ã€Œä¸»æ•°æ®æ–‡ä»¶ã€ã€ŒSMILES CSVã€ã€Œä½¿ç”¨è¯ç‰©ç»“æ„æ¨¡å¼ã€ã€ŒR æ‰§è¡Œæ–¹å¼ã€ã€ŒRscript å‘½ä»¤ã€ã€Œconda.bat è·¯å¾„ã€ã€ŒConda R ç¯å¢ƒåã€ã€Œåˆ†ç»„å­—æ®µï¼ˆgroup_columnï¼‰ã€ã€Œç±»ç¾¤å­—æ®µï¼ˆcluster_columnï¼‰ã€ã€Œéšæœºç§å­ï¼ˆseedï¼‰ã€ã€Œå¾…ç­›é€‰ clustersï¼ˆé€—å·åˆ†éš”ï¼‰ã€ä»¥åŠ gene_sizeã€output_dimã€batch_sizeã€lr æ’°å†™ç®€çŸ­ã€å•ç»†èƒ/LabFlow è¯­å¢ƒä¸‹çš„è¯´æ˜ã€‚
  3. **æ¥å…¥**ï¼šåœ¨ `App.tsx` æ¯ä¸ªå¯¹åº”è¡¨å•é¡¹çš„ label æ–‡å­—åæ’å…¥ `<ParamTooltip text={PARAM_TOOLTIPS["â€¦"]} />`ï¼ˆæˆ– `.gene_size` ç­‰ï¼‰ã€‚
- **æ ·å¼**ï¼š`tokens.css` æ–°å¢ `.param-tooltip-wrap`ã€`.param-tooltip-trigger`ï¼ˆåœ†å½¢ç°åº•ã€Œ?ã€ã€hover æ—¶ accent è‰²ï¼‰ã€`.param-tooltip-bubble`ï¼ˆabsoluteã€æ·±è‰²ã€åœ†è§’ã€ä¸‰è§’ç®­å¤´ï¼‰ã€‚
- **Checkfix**ï¼š`npm run build` é€šè¿‡ï¼›æ— æ–°å¢ lint æŠ¥é”™ã€‚

### [2026-02-11] Seurat è§£æ metadata å€¼æ˜¾ç¤ºï¼šæ·»åŠ ç»†èƒæ•°ç»Ÿè®¡ï¼ˆç±»ä¼¼ R table()ï¼‰
- **ç”¨æˆ·éœ€æ±‚**ï¼šç‚¹å‡» metadata å­—æ®µï¼ˆå¦‚ celltypeï¼‰åï¼Œä¸ä»…è¦æ˜¾ç¤ºæœ‰å“ªäº›åˆ†ç±»ï¼Œè¿˜è¦æ˜¾ç¤ºæ¯ä¸ªåˆ†ç±»æœ‰å¤šå°‘ä¸ªç»†èƒï¼ˆç±»ä¼¼ R çš„ `table()` å‡½æ•°ï¼‰ï¼Œä¾¿äºç”¨æˆ·äº†è§£æ•°æ®åˆ†å¸ƒã€‚
- **é—®é¢˜**ï¼šä¹‹å‰åªè¿”å›å”¯ä¸€å€¼åˆ—è¡¨ï¼Œç”¨æˆ·çœ‹åˆ°"å…± 0 ä¸ª"ï¼Œä¸”æ— æ³•çŸ¥é“æ¯ä¸ªåˆ†ç±»çš„ç»†èƒæ•°ã€‚
- **æ–¹æ¡ˆ**ï¼š
  1. **åç«¯**ï¼ˆ`backend/app/services/seurat_inspector.py`ï¼‰ï¼šä½¿ç”¨ pandas Series çš„ `value_counts()` ç»Ÿè®¡æ¯ä¸ªå€¼çš„è®¡æ•°ï¼ˆç±»ä¼¼ R `table()`ï¼‰ï¼Œè¿”å›ç»“æ„ä» `dict[str, list[str]]` æ”¹ä¸º `dict[str, list[dict[str, Any]]]`ï¼Œæ¯ä¸ªå…ƒç´ åŒ…å« `{"value": str, "count": int}`ï¼›æŒ‰è®¡æ•°é™åºã€å€¼å‡åºæ’åºï¼›è·³è¿‡ NaN/None/ç©ºå­—ç¬¦ä¸²ã€‚
  2. **å‰ç«¯ç±»å‹**ï¼ˆ`frontend/src/services/api.ts`ï¼‰ï¼š`metadata_column_values` ç±»å‹ä» `Record<string, string[]>` æ”¹ä¸º `Record<string, Array<{ value: string; count: number }>>`ã€‚
  3. **å‰ç«¯æ˜¾ç¤º**ï¼ˆ`frontend/src/App.tsx`ï¼‰ï¼šåœ¨ metadata-values-list ä¸­ï¼Œæ¯ä¸ª chip æ˜¾ç¤º `{value} ({count})`ï¼Œcount ç”¨ accent è‰²ã€åŠ ç²—ï¼›hover æ—¶ tooltip æ˜¾ç¤ºå®Œæ•´ä¿¡æ¯ã€‚
- **æ ·å¼**ï¼š`.value-chip` æ”¹ä¸º `inline-flex` æ”¯æŒ gapï¼Œ`.value-count` ç”¨ accent è‰²ã€åŠ ç²—ã€ç¨å°å­—å·ã€‚
- **Checkfix**ï¼šåç«¯ `ruff check` é€šè¿‡ï¼›å‰ç«¯ `npm run build` é€šè¿‡ã€‚

### [2026-02-11] ä¿®å¤ Conda R ç¯å¢ƒåä¸‹æ‹‰é€‰é¡¹ä¸æ˜¾ç¤ºçš„é—®é¢˜
- **ç”¨æˆ·åé¦ˆ**ï¼šConda R ç¯å¢ƒååœ¨å‰ç«¯åˆ·æ–°ä¸å‡ºæ¥é€‰é¡¹ï¼Œåªèƒ½æ‰‹åŠ¨è¾“å…¥ã€‚
- **é—®é¢˜è¯Šæ–­**ï¼š
  1. åˆå§‹åŠ è½½æ—¶è°ƒç”¨ `getCondaEnvs()` æ— å‚æ•°ï¼Œåç«¯å¯èƒ½æ— æ³•è·å–ç¯å¢ƒåˆ—è¡¨ï¼ˆå³ä½¿æ‰¾åˆ°äº† conda.bat candidatesï¼‰
  2. åç«¯é€»è¾‘ï¼šæ—  `conda_bat` å‚æ•°æ—¶ç”¨ç¬¬ä¸€ä¸ª candidate è·å–ç¯å¢ƒåˆ—è¡¨ï¼Œä½†å¯èƒ½æ‰§è¡Œå¤±è´¥è¿”å›ç©ºåˆ—è¡¨
  3. æ‰‹åŠ¨è¾“å…¥ conda.bat è·¯å¾„æ—¶ï¼Œæ²¡æœ‰è‡ªåŠ¨è·å–å¯¹åº”çš„ç¯å¢ƒåˆ—è¡¨
- **ä¿®å¤æ–¹æ¡ˆ**ï¼š
  1. **åˆå§‹åŠ è½½ä¼˜åŒ–**ï¼šè·å–åˆ° `conda_bat_candidates` åï¼Œå¦‚æœ `conda_envs` ä¸ºç©ºï¼Œç”¨ç¬¬ä¸€ä¸ª candidate é‡æ–°è°ƒç”¨ `getCondaEnvs(firstBat)` è·å–ç¯å¢ƒåˆ—è¡¨
  2. **æ‰‹åŠ¨è¾“å…¥å¢å¼º**ï¼šåœ¨æ‰‹åŠ¨è¾“å…¥ conda.bat è·¯å¾„çš„ `onChange` ä¸­ï¼Œå¦‚æœè¾“å…¥ä»¥ "conda.bat" ç»“å°¾ï¼Œè‡ªåŠ¨è°ƒç”¨ `getCondaEnvs()` è·å–ç¯å¢ƒåˆ—è¡¨
- **ä»£ç **ï¼š`frontend/src/App.tsx` ç¬¬ 87-105 è¡Œï¼ˆåˆå§‹åŠ è½½ useEffectï¼‰å’Œç¬¬ 523-535 è¡Œï¼ˆæ‰‹åŠ¨è¾“å…¥ input onChangeï¼‰
- **Checkfix**ï¼šå‰ç«¯ `npm run build` é€šè¿‡ï¼›æ—  lint é”™è¯¯ã€‚

### [2026-02-11] ÑµÁ·ÈÕÖ¾ÂÖÑ¯ 404 ĞŞ¸´ + Ç°¶ËÍ£Ö¹ÑµÁ·ÈÎÎñ
- ÓÃ»§·´À¡£º
  1. ÑµÁ·ÔÚÅÜ£¬µ«Ç°¶Ë³ÖĞøÏÔÊ¾¡°»ñÈ¡ÈÕÖ¾Ê§°Ü¡±¡£
  2. ĞèÒªÔÚÇ°¶ËÖ±½ÓÍ£Ö¹ÑµÁ·ÈÎÎñ¡£
- ¸ùÒò£º
  1. ÈÎÎñÑµÁ·½áÊøÇ°²Å»ØĞ´ `job.log_path`£¬ÂÖÑ¯ `/api/jobs/{id}/log` Ê±³£ÄÃµ½ 404¡£
  2. ºó¶ËÈ±ÉÙÈ¡ÏûÈÎÎñ API£¬Ç°¶ËÒ²Ã»ÓĞ¡°Í£Ö¹ÈÎÎñ¡±Èë¿Ú¡£
- ÊµÏÖ£º
  1. ºó¶Ë `GET /api/jobs/{job_id}/log` ¸ÄÎª¡°ÈÕÖ¾Î´¾ÍĞ÷Ò²·µ»Ø 200 + ¿Õ×Ö·û´®¡±¡£
  2. ºó¶ËĞÂÔö `POST /api/jobs/{job_id}/cancel`£ºqueued Ö±½Ó±ê¼Ç `canceled`£»running ±ê¼Ç `cancel_requested` ²¢³¢ÊÔÖÕÖ¹ÑµÁ·½ø³Ì£¨Windows ÓÃ `taskkill /T /F`£©¡£
  3. ÈÎÎñ¶ÓÁĞÔÚÑµÁ·Æô¶¯Ç°Ô¤Ğ´ `log_path`£¬²¢ÔÚ²¶»ñµ½È¡ÏûÇëÇóÊ±½«×´Ì¬ÊÕÁ²Îª `canceled`¡£
  4. Ç°¶ËĞÂÔö¡°Í£Ö¹ÑµÁ·ÈÎÎñ¡±°´Å¥£¨½ö queued/running ÏÔÊ¾£©£¬½ÓÈë cancel API£»ÈÕÖ¾ÂÖÑ¯Ê§°ÜÊ±²»ÔÙ·´¸´Ğ´¡°»ñÈ¡ÈÕÖ¾Ê§°Ü¡±Õ¼Î»ÎÄ°¸¡£
- ÎÄµµÍ¬²½£º
  1. ĞÂÔö `docs/api/jobs.md`£¨°üº¬ cancel Óë log ĞĞÎª£©¡£
  2. ¸üĞÂ `docs/LabFlowÇ°¶ËÓÃ»§²Ù×÷ËµÃ÷.md`£¨ĞÂÔö¡°Í£Ö¹ÑµÁ·ÈÎÎñ¡±²Ù×÷ËµÃ÷£©¡£
- Checkfix£º
  - ´ıÖ´ĞĞ²¢»ØÌî¡£
- Checkfix ½á¹û»ØÌî£º
  - `python -m black . --check`£ºÊ§°Ü£¬ÌáÊ¾ 2 ¸öÎÄ¼şĞè¸ñÊ½»¯£º`backend/app/api/runtime.py`¡¢`backend/app/services/seurat_inspector.py`¡£
  - `npm run lint`£¨frontend£©£ºÍ¨¹ı¡£
  - `npm run build`£¨frontend£©£ºÍ¨¹ı¡£
  - ºó¶Ë pytest£ºµ±Ç°»·¾³È±ÉÙ pytest£¨`No module named pytest`£©£»`uv run pytest` ÊÜ±¾»úÈ¨ÏŞÏŞÖÆ£¨¾Ü¾ø·ÃÎÊ£©¡£
- [2026-02-11] Checkfix ×·¼Ó£¨code-debugger£©
  - ÓÃ»§±¾»úÈ·ÈÏ£º`python -m black backend/app/services/seurat_inspector.py --check` Í¨¹ı¡£
  - ÒÑÍê³É×Ô¶¯¸ñÊ½»¯ÊÕÁ²£º`backend/app/api/runtime.py`¡¢`backend/app/api/jobs.py`£¨Í³Ò»¸ñÊ½Óë»»ĞĞ£©¡£
  - ¸´¼ìÍ¨¹ı£º
    - `python -m black backend/app/api/runtime.py backend/app/services/seurat_inspector.py --check`
    - `ruff format --check backend/app backend/tests`
    - `ruff check backend/app backend/tests`

### [2026-02-11] ¶¯»­Ê×Ò³ + ¿ªÊ¼·ÖÎöÈë¿Ú£¨¿ÆÑĞÊÓ¾õÉıÎ¬£©
- ÒÕÊõÖ¸µ¼£º
  - Mood: ÀíĞÔ¡¢¾«ÃÜ¡¢¿É¹Û²â¡£
  - Metaphor: ¡°ÊµÑéÌ¨ÉÏµÄ¹âÑ§É¨ÃèÆÁ¡±£¬ÏÈÕ¹Ê¾ÏµÍ³ÄÜÁ¦£¬ÔÙ½øÈë·ÖÎöÁ÷³Ì¡£
- ÊÓ¾õÉó¼ÆÓë²ßÂÔ£º
  1. ¿Õ¼ä£ºÔ­Ò³Ãæ´ò¿ª¼´½øÈë¶àÃæ°å£¬ĞÅÏ¢ÃÜ¶È¸ß£»ĞÂÔö Landing Ê×ÆÁ£¬·Ö²ãÕ¹Ê¾ÄÜÁ¦µãÓëÈë¿Ú¶¯×÷¡£
  2. ÕÅÁ¦£º¼ÓÈëÉ¨Ãè¸ß¹â¡¢Íø¸ñÎÆÀí¡¢¿¨Æ¬ºôÎü¶¯»­£¬ÔöÇ¿¿ÆÑĞÒÇ±í¸Ğ¡£
  3. ÖÊ¸Ğ£º¶à²ã½¥±ä + °ëÍ¸Ã÷¿¨Æ¬ + ÇáÒõÓ°£¬±ÜÃâ´¿Æ½½çÃæ¡£
  4. Î¢½»»¥£ºÈë¿Ú°´Å¥ hover ÌáÉı¡¢Á÷³ÌÒ³½øÈëÊ±ÉÏ¸¡¹ı¶É¡£
- ÊµÊ©¼ÇÂ¼£º
  - `frontend/src/App.tsx`£ºĞÂÔö `hasEntered` Èë¿Ú×´Ì¬£»¼ÓÈë¶¯»­Ê×Ò³ DOM£»µã»÷¡°¿ªÊ¼·ÖÎö¡±ºóäÖÈ¾Ô­ÓĞ 1~7 Á÷³Ì¡£
  - `frontend/src/styles/tokens.css`£ºĞÂÔö landing Ïà¹Ø token ÓëÑùÊ½£¨Íø¸ñ²ã¡¢É¨Ãè¶¯»­¡¢¿¨Æ¬ºôÎü¡¢Èë¿Ú°´Å¥¶¯Ğ§¡¢ÒÆ¶¯¶ËÊÊÅä£©¡£
- ÓÃ»§ËµÃ÷Êé¸üĞÂ£º
  - `docs/LabFlowÇ°¶ËÓÃ»§²Ù×÷ËµÃ÷.md` ĞÂÔö¡°¶¯»­Ê×Ò³Èë¿Ú£¨¿ªÊ¼·ÖÎö£©¡±ÕÂ½Ú£¨Ä¿µÄ¡¢²Ù×÷¡¢Ô¤ÆÚ¡¢³£¼ûÎÊÌâ¡¢»Ø¹ö£©¡£
- ²¿ÊğÎÄµµÁª¶¯¼ì²é£º
  - ÒÑ¼ì²é `docs/²¿ÊğÎÄµµ.md`£¬±¾´Î½öÇ°¶ËÊÓ¾õÓëÈë¿Ú½»»¥±ä»¯£¬²»Éæ¼°²¿ÊğÃüÁîÓë»·¾³±äÁ¿£¬ÎŞĞè¸Ä¶¯¡£
- Checkfix£º
  - `npm run lint`£¨frontend£©Í¨¹ı¡£
  - `npm run build`£¨frontend£©Í¨¹ı¡£

### [2026-02-11] Ê×Ò³×¢²á/µÇÂ¼ + ÓÃ»§ËµÃ÷ÊéÈë¿Ú£¨ÄÚÍøÇáÁ¿ÈÏÖ¤£©
- Ä¿±ê£ºÔÚ¶¯»­Ê×Ò³¼ÓÈë¿ÉÖ±½ÓÊ¹ÓÃµÄÕËºÅÈë¿Ú£¬µÇÂ¼ºó½øÈë·ÖÎöÁ÷³Ì£¬²¢°ÑÓÃ»§ËµÃ÷Êé°´Å¥·ÅÔÚÊ×Ò³¡£
- ºó¶ËÊµÏÖ£¨API-First£©£º
  - ĞÂÔö `backend/app/services/auth_service.py`£ºSQLite ±¾µØ¿â + PBKDF2-SHA256 ÃÜÂë¹şÏ£ + session token¡£
  - ĞÂÔö `backend/app/api/auth.py`£º`/api/auth/register`¡¢`/api/auth/login`¡¢`/api/auth/me`¡¢`/api/auth/logout`¡¢`/api/auth/user-guide`¡£
  - ĞÂÔö `backend/app/auth.py`£ºBearer ½âÎöÓëÈÏÖ¤ÒÀÀµ¡£
  - `backend/app/main.py` Â·ÓÉ¹ÒÔØ£ºÒµÎñ API Í³Ò»×ß `require_auth`£¨¿ÉÓÉ `LABFLOW_AUTH_REQUIRED` ¿ØÖÆ£©¡£
  - `backend/app/core/config.py` Ôö¼ÓÈÏÖ¤ÅäÖÃ£º`LABFLOW_AUTH_REQUIRED`¡¢`LABFLOW_AUTH_DB_PATH`¡¢`LABFLOW_AUTH_SESSION_TTL_HOURS`¡£
- Ç°¶ËÊµÏÖ£¨Ê×Ò³Èë¿ÚÓë½»»¥£©£º
  - `frontend/src/App.tsx`£ºÊ×Ò³ĞÂÔö×¢²á/µÇÂ¼±íµ¥¡¢µÇÂ¼Ì¬ÏÔÊ¾¡¢ÍË³öµÇÂ¼¡¢Î´µÇÂ¼½ûÓÃ¡°¿ªÊ¼·ÖÎö¡±¡£
  - Ê×Ò³¡°ÓÃ»§ËµÃ÷Êé¡±°´Å¥¸ÄÎªºó¶ËÔÚÏßÎÄµµ½Ó¿Ú£º`/api/auth/user-guide`¡£
  - `frontend/src/services/api.ts`£ºĞÂÔö auth API ·â×°Óë token ±¾µØ´æÈ¡£»ÇëÇó×Ô¶¯Ğ¯´ø Bearer Token¡£
  - `frontend/src/styles/tokens.css`£ºĞÂÔö auth panel¡¢tabs¡¢°´Å¥ÑùÊ½¡£
- ÎÄµµÍ¬²½£º
  - ĞÂÔö `docs/api/auth.md`¡£
  - ¸üĞÂ `docs/LabFlowÇ°¶ËÓÃ»§²Ù×÷ËµÃ÷.md`£¨ĞÂÔöÊ×Ò³×¢²á/µÇÂ¼²½ÖèÓë¹æÔò£©¡£
  - ¸üĞÂ `docs/²¿ÊğÎÄµµ.md`£¨ĞÂÔöÈÏÖ¤Ïà¹Ø»·¾³±äÁ¿£©¡£
  - ¸üĞÂ `README.md`£¨ĞÂÔöµÇÂ¼ÈÏÖ¤ËµÃ÷Èë¿Ú£©¡£
- ÓÃ»§¿ÉÀí½âĞÔ¼ì²é½áÂÛ£º
  - ½áÂÛ£º»ù±¾¿ÉÈÃĞÂÊÖÉÏÊÖ¡£ËµÃ÷ÊéÒÑ¸²¸Ç¡°ÏÈµÇÂ¼ÔÙ¿ªÊ¼·ÖÎö¡±µÄÂ·¾¶¡¢ÕËºÅ¹æÔò¡¢³£¼û±¨´í¡£
  - ÈÔ½¨Òé£¨ºóĞø¿É×ö£©£º²¹Ò»ÕÅ¡°×¢²á->µÇÂ¼->¿ªÊ¼·ÖÎö¡±Á÷³Ì½ØÍ¼£¬½øÒ»²½½µµÍÊ×´ÎÈÏÖª³É±¾¡£
- Checkfix£º
  - `ruff check backend/app backend/tests` Í¨¹ı¡£
  - `ruff format --check backend/app backend/tests` Í¨¹ı¡£
  - `npm run lint`£¨frontend£©Í¨¹ı¡£
  - `npm run build`£¨frontend£©Í¨¹ı¡£
  - `python -m pytest -q backend/tests/test_auth_api.py` Î´Ö´ĞĞ£¨µ±Ç°»·¾³È±ÉÙ pytest£©¡£
  - ÔËĞĞÊ± smoke£¨ÉèÖÃ `LABFLOW_AUTH_DB_PATH=.debug/auth_smoke.db`£©Í¨¹ı£º
    - `POST /api/auth/register` -> 200
    - `GET /api/auth/me` -> 200
    - `POST /api/auth/logout` -> 200
    - logout ºó `GET /api/auth/me` -> 401
    - `GET /api/auth/user-guide` -> 200
    - Î´µÇÂ¼·ÃÎÊ `GET /api/jobs` -> 401

### [2026-02-11 16:31 +08:00] ÈÎÎñÂÖÑ¯ÌåÑéÔöÇ¿£ºGPU ÊµÊ±¿´°å + µÇÂ¼Ö±´ïÁ÷³Ì + ÈÕÖ¾ÔçÏÔ
- ÎÊÌâÃèÊö
  - ÓÃ»§Ï£ÍûÂÖÑ¯½×¶ÎÓĞ¸üÖ±¹Û·´À¡£¨¶¯»­»òÏÔ¿¨×´Ì¬£©£¬²¢ÇÒÈÕÖ¾²»ÒªµÈ¡°Í£Ö¹ÈÎÎñ¡±ºó²Å³öÏÖ¡£
  - µÇÂ¼ºóÈÔĞèÔÙµã¡°¿ªÊ¼·ÖÎö¡±£¬Á÷³Ì²»¹»Ë³»¬¡£
  - Êı¾İÌá½»Ò³È±ÉÙÓÃ»§ËµÃ÷ÊéÈë¿Ú¡£
- ¸ùÒò¶¨Î»
  - Ç°¶ËÎ´½ÓÈë GPU ÂÖÑ¯Êı¾İÔ´¡£
  - µÇÂ¼³É¹¦Ö»ÉèÖÃ token/user£¬Î´ÇĞ»» `hasEntered=true`¡£
  - ÈÎÎñÈÕÖ¾½Ó¿ÚÓÅÏÈ¶ÁÈ¡ `train.log`£¬ÑµÁ·ÔçÆÚ¸ÃÎÄ¼ş¿ÉÄÜÎª¿Õ£¬µ¼ÖÂÇ°¶ËÏÔÊ¾¡°µÈ´ıÈÕÖ¾¡±¡£
- ½â¾ö·½°¸
  - ºó¶ËĞÂÔö `GET /api/runtime/gpu-stats`£¨`nvidia-smi` ½âÎö£©¡£
  - Ç°¶ËÔÚÈÎÎñ `running` Ê±ÂÖÑ¯ GPU Ö¸±ê²¢Õ¹Ê¾ `GPU Runtime` ¿¨Æ¬£¨º¬Âö³å¶¯»­ºÍ×ÊÔ´Ìõ£©¡£
  - µÇÂ¼³É¹¦ºó×Ô¶¯ `setHasEntered(true)` Ö±´ï·ÖÎöÒ³¡£
  - ÔÚ `1) ÉÏ´«Êı¾İ` Ãæ°åĞÂÔö `User Guide` °´Å¥¡£
  - `GET /api/jobs/{job_id}/log` Ôö¼Ó fallback£ºµ± `job.log_path` ÄÚÈİÎª¿ÕÊ±£¬»ØÍË¶ÁÈ¡ `backend/artifacts/jobs/{job_id}/logger/log.txt`¡£
  - Í£Ö¹ÑµÁ·°´Å¥½öÔÚ `running` ÏÔÊ¾£¬²»ÔÚ `queued` ÏÔÊ¾¡£
- ´úÂë±ä¸ü£¨ÎÄ¼ş£©
  - `backend/app/api/runtime.py`
  - `backend/app/api/jobs.py`
  - `frontend/src/services/api.ts`
  - `frontend/src/App.tsx`
  - `frontend/src/styles/tokens.css`
  - `docs/api/runtime.md`£¨ĞÂÔö£©
  - `docs/api/auth.md`
  - `docs/api/jobs.md`
  - `docs/LabFlowÇ°¶ËÓÃ»§²Ù×÷ËµÃ÷.md`
- Checkfix ½á¹û
  - `ruff format --check backend/app backend/tests` -> Í¨¹ı£¨ÏÈ×Ô¶¯¸ñÊ½»¯ `backend/app/api/auth.py`£©
  - `ruff check backend/app backend/tests` -> Í¨¹ı
  - `npm run lint`£¨frontend£©-> Í¨¹ı
  - `npm run build`£¨frontend£©-> Í¨¹ı
- ²¿Êğ/ÎÄµµÁª¶¯¼ì²é
  - ÒÑ¼ì²é²¿ÊğÎÄµµ£º±¾´Î²»ĞÂÔö²¿ÊğÃüÁî£¬½öĞÂÔöÔËĞĞÊ± API ÓëÇ°¶ËÕ¹Ê¾£¬ÎŞĞèĞŞ¸Ä `docs/²¿ÊğÎÄµµ.md`¡£

### [2026-02-11 16:45 +08:00] ä¿®å¤ï¼š/api/auth/user-guide æ‰“ä¸å¼€ + æ–‡æ¡£ä¹±ç  + å‰ç«¯è¿æ¥æ‹’ç»
- é—®é¢˜ç°è±¡
  - å‰ç«¯æ§åˆ¶å°å¤§é‡ `ERR_CONNECTION_REFUSED`ï¼ˆ`/api/health`ã€`/api/auth/me`ã€`/api/auth/register`ï¼‰ã€‚
  - `http://localhost:8000/api/auth/user-guide` æ— æ³•è®¿é—®ã€‚
  - `docs/LabFlowå‰ç«¯ç”¨æˆ·æ“ä½œè¯´æ˜.md` ååŠæ®µä¹±ç ã€‚
- æ ¹å› 
  1. åç«¯å¯åŠ¨å¤±è´¥ï¼š`backend/app/api/auth.py` ä¸­ `user_guide` è¿”å›ç±»å‹æ ‡æ³¨ä¸º `FileResponse | HTMLResponse`ï¼ŒFastAPI åœ¨è·¯ç”±å»ºæ¨¡æ—¶æŠ¥é”™å¹¶é€€å‡ºã€‚
  2. ç”¨æˆ·è¯´æ˜ä¹¦æ–‡ä»¶ç¼–ç å·²æŸåï¼Œå¯¼è‡´ `/api/auth/user-guide` è¿”å› `User guide encoding is not supported`ã€‚
- ä¿®å¤
  1. `backend/app/api/auth.py`
     - `@router.get("/user-guide", response_model=None)`
     - è¿”å›ç±»å‹æ”¹ä¸º `Response`ï¼Œé¿å… FastAPI å°†å“åº”ç±» Union å½“æˆ Pydantic å­—æ®µã€‚
  2. å½»åº•é‡å»º `docs/LabFlowå‰ç«¯ç”¨æˆ·æ“ä½œè¯´æ˜.md`
     - ä»¥ UTF-8 é‡å†™æ•´ä»½æ–‡æ¡£ï¼Œè¦†ç›–ç™»å½•ã€ä¸Šä¼ ã€æ ¡éªŒã€è®­ç»ƒã€GPU çœ‹æ¿ã€å¸¸è§æ•…éšœä¸å›æ»šã€‚
- éªŒè¯
  - æœ¬åœ°å¯åŠ¨åç«¯å¹¶è°ƒç”¨ï¼š
    - `GET /api/health` -> 200
    - `GET /api/auth/user-guide` -> 200ï¼ˆHTMLï¼‰
    - `GET /api/auth/user-guide?raw=true` -> 200ï¼ˆMarkdownï¼‰
  - Checkfixï¼š
    - `ruff format --check backend/app backend/tests` é€šè¿‡
    - `ruff check backend/app backend/tests` é€šè¿‡
    - `npm run lint` é€šè¿‡
    - `npm run build` é€šè¿‡
- ç”¨æˆ·ä¾§æ“ä½œæç¤º
  - å‡ºç° `ERR_CONNECTION_REFUSED` æ—¶å…ˆç¡®è®¤åç«¯å·²è¿è¡Œåœ¨ `8000` ç«¯å£ã€‚
