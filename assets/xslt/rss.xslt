<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" >
<xsl:output method="html" encoding="utf-8" />
<xsl:template match="/rss">
	<xsl:text disable-output-escaping="yes">&lt;!DOCTYPE html &gt;</xsl:text>
	<html>
	<head>
		<xsl:text disable-output-escaping="yes"><![CDATA[
		<meta charset="utf-8">
	<meta name="viewport" content="width=device-width, initial-scale=1.0">
	<title>RSS Feed (Styled)</title>
	<link rel="stylesheet" type="text/css" href="http://redbkit.github.io/redbKIT/assets/css/styles_feeling_responsive.css">
	<script src="http://redbkit.github.io/redbKIT/assets/js/modernizr.min.js"></script>

	<script src="https://ajax.googleapis.com/ajax/libs/webfont/1.5.18/webfont.js"></script>
	<script>
		WebFont.load({
			google: {
				families: [ 'Lato:400,700,400italic:latin', 'Volkhov::latin' ]
			}
		});
	</script>

	<noscript>
		<link href='http://fonts.googleapis.com/css?family=Lato:400,700,400italic%7CVolkhov' rel='stylesheet' type='text/css'>
	</noscript>


	<!-- Search Engine Optimization -->
	<meta name="description" content="»redbKIT:« a MATLAB library for reduced-order modeling of parametrized PDEs">
  	
	
	
	
	


	<!-- Facebook Open Graph -->
	<meta property="og:title" content="RSS Feed (Styled)">
	<meta property="og:description" content="»redbKIT:« a MATLAB library for reduced-order modeling of parametrized PDEs">
	<meta property="og:url" content="http://redbkit.github.io/redbKIT/assets/xslt/rss.xslt">
	<meta property="og:locale" content="en_EN">
	<meta property="og:type" content="website">
	<meta property="og:site_name" content="redbKIT">
	
	


	

	<link type="text/plain" rel="author" href="http://redbkit.github.io/redbKIT/humans.txt">

	

	

	<link rel="icon" sizes="32x32" href="http://redbkit.github.io/redbKIT/assets/img/favicon-32x32.png">

	<link rel="icon" sizes="192x192" href="http://redbkit.github.io/redbKIT/assets/img/touch-icon-192x192.png">

	<link rel="apple-touch-icon-precomposed" sizes="180x180" href="http://redbkit.github.io/redbKIT/assets/img/apple-touch-icon-180x180-precomposed.png">

	<link rel="apple-touch-icon-precomposed" sizes="152x152" href="http://redbkit.github.io/redbKIT/assets/img/apple-touch-icon-152x152-precomposed.png">

	<link rel="apple-touch-icon-precomposed" sizes="144x144" href="http://redbkit.github.io/redbKIT/assets/img/apple-touch-icon-144x144-precomposed.png">

	<link rel="apple-touch-icon-precomposed" sizes="120x120" href="http://redbkit.github.io/redbKIT/assets/img/apple-touch-icon-120x120-precomposed.png">

	<link rel="apple-touch-icon-precomposed" sizes="114x114" href="http://redbkit.github.io/redbKIT/assets/img/apple-touch-icon-114x114-precomposed.png">

	
	<link rel="apple-touch-icon-precomposed" sizes="76x76" href="http://redbkit.github.io/redbKIT/assets/img/apple-touch-icon-76x76-precomposed.png">

	<link rel="apple-touch-icon-precomposed" sizes="72x72" href="http://redbkit.github.io/redbKIT/assets/img/apple-touch-icon-72x72-precomposed.png">

	<link rel="apple-touch-icon-precomposed" href="http://redbkit.github.io/redbKIT/assets/img/apple-touch-icon-precomposed.png">	

	<meta name="msapplication-TileImage" content="http://redbkit.github.io/redbKIT/assets/img/msapplication_tileimage.png">

	<meta name="msapplication-TileColor" content="#fabb00">


	

		]]></xsl:text>
	</head>
	<body id="top-of-page">
		<xsl:text disable-output-escaping="yes"><![CDATA[
		<div id="navigation" class="sticky">
  <nav class="top-bar" role="navigation" data-topbar>
    <ul class="title-area">
      <li class="name">
      <h1 class="show-for-small-only"><a href="http://redbkit.github.io/redbKIT" class="icon-tree"> redbKIT</a></h1>
    </li>
       <!-- Remove the class "menu-icon" to get rid of menu icon. Take out "Menu" to just have icon alone -->
      <li class="toggle-topbar menu-icon"><a href="#"><span>Navigation</span></a></li>
    </ul>
    <section class="top-bar-section">

      <ul class="right">
        

              

          
          
        

              

          
          
        

              

          
          
        

              

          
          
        

              

          
          
        

              

          
          
        

              

          
          
            
            
              <li class="divider"></li>
              <li><a href="http://redbkit.github.io/redbKIT/search/">Search</a></li>

            
            
          
        

              

          
          
            
            
              <li class="divider"></li>
              <li><a href="http://redbkit.github.io/redbKIT/contact/">Contact</a></li>

            
            
          
        
        
      </ul>

      <ul class="left">
        

              

          
          

            
            
              <li><a href="http://redbkit.github.io/redbKIT/">About</a></li>
              <li class="divider"></li>

            
            
          
        

              

          
          

            
            
              <li><a href="http://redbkit.github.io/redbKIT/download-install/">Download and Install</a></li>
              <li class="divider"></li>

            
            
          
        

              

          
          

            
            
              <li><a href="http://redbkit.github.io/redbKIT/getting-started/">Getting Started</a></li>
              <li class="divider"></li>

            
            
          
        

              

          
          

            
            
              <li><a href="http://redbkit.github.io/redbKIT/probgallery/">Problems Gallery</a></li>
              <li class="divider"></li>

            
            
          
        

              

          
          

            
            
              <li><a href="http://redbkit.github.io/redbKIT/contributors/">Development</a></li>
              <li class="divider"></li>

            
            
          
        

              

          
          

            
            

              <li class="has-dropdown">
                <a href="http://redbkit.github.io/redbKIT/documentation/">Documentation</a>

                  <ul class="dropdown">
                    

                      


			
                      
		      
			<li class="has-dropdown">
                <a href="http://redbkit.github.io/redbKIT/math/ADR/">Advection diffusion reaction</a>

                  <ul class="dropdown">
                    

                      

                      <li><a href="http://redbkit.github.io/redbKIT/math/AdvDiffReact/">Advection diffusion reaction equations</a></li>
			

                      

                      <li><a href="http://redbkit.github.io/redbKIT/math/AdvDiffReactEx/">Simulation Setup</a></li>
			
</ul>
 </li>
		                          
		      
			

                      


			
                      
		      
			<li class="has-dropdown">
                <a href="http://redbkit.github.io/redbKIT/math/CFD/">Computational Fluid Dynamics</a>

                  <ul class="dropdown">
                    

                      

                      <li><a href="http://redbkit.github.io/redbKIT/math/NavierStokes/">Navier-Stokes equations</a></li>
			

                      

                      <li><a href="http://redbkit.github.io/redbKIT/math/NavierStokesEx/">Simulation Setup</a></li>
			
</ul>
 </li>
		                          
		      
			

                      


			
                      
		      
			<li class="has-dropdown">
                <a href="http://redbkit.github.io/redbKIT/math/CSM/">Computational Solid Mechanics</a>

                  <ul class="dropdown">
                    

                      

                      <li><a href="http://redbkit.github.io/redbKIT/math/SolidMechanics/">Hyperelasticity: static and dynamics</a></li>
			

                      

                      <li><a href="http://redbkit.github.io/redbKIT/math/SolidMechanicsEx/">Simulation Setup</a></li>
			
</ul>
 </li>
		                          
		      
			

                      


			
                      
		      
			<li class="has-dropdown">
                <a href="http://redbkit.github.io/redbKIT/math/FSI/">Fluid-Structure Interaction</a>

                  <ul class="dropdown">
                    

                      

                      <li><a href="http://redbkit.github.io/redbKIT/math/FSIsolver/">Monolithic FSI solver</a></li>
			

                      

                      <li><a href="http://redbkit.github.io/redbKIT/math/FSIEx/">Simulation Setup</a></li>
			

                      

                      <li><a href="http://redbkit.github.io/redbKIT/math/FSIprestress/">Tissue prestress</a></li>
			
</ul>
 </li>
		                          
		      
			

                      


			
                      
		      
                      <li><a href="http://redbkit.github.io/redbKIT/math/essentialbc/">Essential Boundary Conditions</a></li>
		                          
		      
			
                  </ul>

              </li>
              <li class="divider"></li>
            
          
        

              

          
          
        

              

          
          
        
        
      </ul>
    </section>
  </nav>
</div><!-- /#navigation -->

		



<div id="masthead-with-text">
	<div class="row">
		<div class="medium-9 small-centered columns">
			<a id="forkme_banner" href="https://github.com/redbKIT/redbKIT">View on GitHub</a>
			<section id="downloads">
             			 <a class="zip_download_link" href="https://github.com/redbKIT/redbKIT/archive/v2.1.zip">Download this project as a .zip file</a>
             			 <a class="tar_download_link" href="https://github.com/redbKIT/redbKIT/archive/v2.1.tar.gz">Download this project as a tar.gz file</a>
			</section>
			<div class="masthead-title"><b>redbKIT</b> <br/> <small>a MATLAB library for reduced-order modeling of PDEs</small></div>
		</div> 
	</div> 
</div> 







		


<div class="alert-box warning text-center"><p>This <a href="https://en.wikipedia.org/wiki/RSS" target="_blank">RSS feed</a> is meant to be used by <a href="https://en.wikipedia.org/wiki/Template:Aggregators" target="_blank">RSS reader applications and websites</a>.</p>
</div>



		]]></xsl:text>
		<header class="t30 row">
	<p class="subheadline"><xsl:value-of select="channel/description" disable-output-escaping="yes" /></p>
	<h1>
		<xsl:element name="a">
			<xsl:attribute name="href">
				<xsl:value-of select="channel/link" />
			</xsl:attribute>
			<xsl:value-of select="channel/title" disable-output-escaping="yes" />
		</xsl:element>
	</h1>
</header>
<ul class="accordion row" data-accordion="">
	<xsl:for-each select="channel/item">
		<li class="accordion-navigation">
			<xsl:variable name="slug-id">
				<xsl:call-template name="slugify">
					<xsl:with-param name="text" select="guid" />
				</xsl:call-template>
			</xsl:variable>
			<xsl:element name="a">
				<xsl:attribute name="href"><xsl:value-of select="concat('#', $slug-id)"/></xsl:attribute>
				<xsl:value-of select="title"/>
				<br/>
				<small><xsl:value-of select="pubDate"/></small>
			</xsl:element>
			<xsl:element name="div">
				<xsl:attribute name="id"><xsl:value-of select="$slug-id"/></xsl:attribute>
				<xsl:attribute name="class">content</xsl:attribute>
				<h1>
					<xsl:element name="a">
						<xsl:attribute name="href"><xsl:value-of select="link"/></xsl:attribute>
						<xsl:value-of select="title"/>
					</xsl:element>
				</h1>
				<xsl:value-of select="description" disable-output-escaping="yes" />
			</xsl:element>
		</li>
	</xsl:for-each>
</ul>

		<xsl:text disable-output-escaping="yes"><![CDATA[
		    <div id="up-to-top" class="row">
      <div class="small-12 columns" style="text-align: right;">
        <a class="iconfont" href="#top-of-page">&#xf108;</a>
      </div><!-- /.small-12.columns -->
    </div><!-- /.row -->

      <div id="subfooter">
        <nav class="row">
          <section id="subfooter-left" class="small-12 medium-6 columns credits">
            <p><a href="http://cmcs.epfl.ch/people/negri">&copy;Federico Negri</a>. Created with <a href="http://jekyllrb.com/" target="_blank">Jekyll</a> based on <a href="http://phlow.github.io/feeling-responsive/">Feeling Responsive</a>.</p>
          </section>

          <section id="subfooter-right" class="small-12 medium-6 columns">
            <ul class="inline-list social-icons">
            
            </ul>
          </section>
        </nav>
      </div><!-- /#subfooter -->
    </footer>

		


<script src="http://redbkit.github.io/redbKIT/assets/js/javascript.min.js"></script>







<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-64466050-1', 'auto');
  ga('set', 'anonymizeIp', true);
  ga('send', 'pageview');

</script>








		]]></xsl:text>
	</body>
	</html>
</xsl:template>
<xsl:template name="slugify">
	<xsl:param name="text" select="''" />
	<xsl:variable name="dodgyChars" select="' ,.#_-!?*:;=+|&amp;/\\'" />
	<xsl:variable name="replacementChar" select="'-----------------'" />
	<xsl:variable name="lowercase" select="'abcdefghijklmnopqrstuvwxyz'" />
	<xsl:variable name="uppercase" select="'ABCDEFGHIJKLMNOPQRSTUVWXYZ'" />
	<xsl:variable name="lowercased"><xsl:value-of select="translate( $text, $uppercase, $lowercase )" /></xsl:variable>
	<xsl:variable name="escaped"><xsl:value-of select="translate( $lowercased, $dodgyChars, $replacementChar )" /></xsl:variable>
	<xsl:value-of select="$escaped" />
</xsl:template>
</xsl:stylesheet>
