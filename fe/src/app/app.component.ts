import {Component, OnInit} from '@angular/core';
import {MatSidenavModule} from "@angular/material/sidenav";
import {MatToolbarModule} from "@angular/material/toolbar";
import {MatButtonModule} from "@angular/material/button";
import {MatIconModule} from "@angular/material/icon";
import {RouterLink, RouterOutlet} from "@angular/router";
import {BreakpointObserver, Breakpoints} from "@angular/cdk/layout";
import {MatListItem, MatNavList} from "@angular/material/list";

@Component({
  selector: 'app-root',
  imports: [
    MatSidenavModule,
    MatToolbarModule,
    MatButtonModule,
    MatIconModule,
    RouterLink,
    RouterOutlet,
    MatListItem,
    MatNavList,
  ],
  templateUrl: './app.component.html',
  styleUrl: './app.component.scss'
})
export class AppComponent implements OnInit {
  showSidebar = false;
  footerClass = "footer_normal"
  constructor(private responsive: BreakpointObserver) {
  }

  ngOnInit() {
    this.responsive.observe([
      Breakpoints.TabletPortrait,
      Breakpoints.TabletLandscape,
      Breakpoints.HandsetPortrait,
      Breakpoints.HandsetLandscape
    ])
    .subscribe(result => {
      this.showSidebar = false;
      this.footerClass = "footer_normal";
      const breakpoints = result.breakpoints;
      if (
        breakpoints[Breakpoints.TabletPortrait] ||
        breakpoints[Breakpoints.TabletLandscape] ||
        breakpoints[Breakpoints.HandsetPortrait] ||
        breakpoints[Breakpoints.HandsetLandscape]
      ) {
        this.showSidebar = true;
        this.footerClass = "footer_small";
      }
    });
  }
}
