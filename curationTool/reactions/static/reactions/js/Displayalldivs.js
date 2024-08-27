function displayDivs(reactionData) {
    loadAtomMappingDiv(reactionData);
    loadChemInfoDiv(reactionData);
    loadMetaboliteInfoDiv(reactionData);
    refreshSideButtons();
    updateStatusDots('substratesDiv', reactionData.subs_found, reactionData.subs_miriams);
    updateStatusDots('productsDiv', reactionData.prod_found, reactionData.prod_miriams);
    displayReactionMessage(reactionData);
    displayreactioninfo(reactionData);
}
