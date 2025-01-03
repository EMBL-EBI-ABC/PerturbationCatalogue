import {Component, OnInit} from '@angular/core';
import {MaveDBData} from "../../model/mavedb";
import {ActivatedRoute} from "@angular/router";
import {MatCardModule} from "@angular/material/card";
import {MatDividerModule} from "@angular/material/divider";

@Component({
  selector: 'app-details',
  imports: [
    MatCardModule,
    MatDividerModule,
  ],
  templateUrl: './details.component.html',
  styleUrl: './details.component.scss'
})
export class DetailsComponent implements OnInit {
  maveDBRecord!: MaveDBData;

  constructor(private route: ActivatedRoute) { }

  ngOnInit() {
    this.maveDBRecord = this.route.snapshot.data['record']['results'][0];
  }
}
