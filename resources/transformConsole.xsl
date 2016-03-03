<?xml version="1.0" encoding="ISO-8859-1"?>

<xsl:stylesheet version = '1.0'
     xmlns:xsl='http://www.w3.org/1999/XSL/Transform'>
     
	<xsl:output method="text"/>
	
	<xsl:template match="/">
		<xsl:apply-templates select="//module/parameters"/>

		<!-- xsl:apply-templates select="//module/parameters/parCat/comments"/ -->
    </xsl:template>  
    	
	<xsl:template match="module/parameters">
		<!--xsl:apply-templates select="parCat"/>
		<xsl:text>
</xsl:text -->
		<xsl:apply-templates select=".//par"/>
		
		<xsl:text>
</xsl:text>
	</xsl:template>
	
	<xsl:template match="parCat">
		<xsl:apply-templates select="par"/>
	</xsl:template>

	
	<xsl:template match="par">
		<xsl:apply-templates select="par[@console]"/>
		<xsl:apply-templates select="par[@system]"/>
		<!-- xsl:value-of select="normalize-space()"/>	
		<xsl:text> 
</xsl:text -->
	
	</xsl:template>
	
	<xsl:template match="par[@console]">
		<xsl:value-of select="normalize-space()"/>	
		<xsl:text> 
</xsl:text>
	</xsl:template>
	
	<!-- xsl:template match="par[last()]">
		<xsl:value-of select="normalize-space()"/>	
		<xsl:text> </xsl:text>
		<xsl:apply-templates select="../comments"/>
		<xsl:text>
</xsl:text>
		
	</xsl:template -->
		
	<xsl:template match="par[@type='cdata']">
	<!-- ýëåìåíòû cdata âûâîäÿòñÿ áåç íîðìàëèçàöèè-->
		<xsl:value-of select="text()"/>
	</xsl:template>
	    
    
   <!-- <xsl:template match="/worker/module/parameters">
		<xsl:apply-templates select="./par[2]"/>
    </xsl:template> -->

</xsl:stylesheet>

